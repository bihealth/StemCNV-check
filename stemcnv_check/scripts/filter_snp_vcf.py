import vcfpy
from loguru import logger
#from snakemake.script import snakemake

logger.remove()
logger.add(snakemake.log['err'], level='DEBUG')
logger.info(f"Reading in '{snakemake.input[0]}' and will write to '{snakemake.output[0]}'")

def get_filter_settings(filtersetname, config):

    filterset = vcfpy.OrderedDict()
    config_filterset = config['settings']['probe-filter-sets'].get(filtersetname)
    if not config_filterset:
        raise ValueError(f"Filter set {filtersetname} not found in config file.")

    GT_min = config_filterset.get('GenTrainScore')
    if GT_min:
        filterset['GT'] = (f'GT{GT_min}', GT_min, f'GenTrain score below {GT_min}')
    GC_min = config_filterset.get('GenCallScore')
    if GC_min:
        filterset['GC'] = (f'IGC{GC_min}', GC_min, f'GenCall score (IGC) below {GC_min}')

    nonuniq = config_filterset.get('Position.duplicates')
    if nonuniq and nonuniq != 'keep':
        desc = 'Remove records with non-unique positions'
        if nonuniq != 'remove':
            desc += f', but keep one per position with {nonuniq}'
        filterset['nonuniq'] = ('NON-UNIQ', nonuniq, desc)

    return filterset


#TODO: fail on multi-sample vcf or [current] filter if ANY samples is below thresholds IGC ?
def filter_snp_vcf(filterset, input='/dev/stdin', output='/dev/stdout'):

    reader = vcfpy.Reader.from_path(input)

    for _, filtercond in filterset.items():
        reader.header.add_filter_line(
            vcfpy.OrderedDict([
                ('ID', filtercond[0]),
                ('Description', filtercond[2])
            ])
        )

    writer = vcfpy.Writer.from_path(output, reader.header)

    def apply_uniq_pos_filter(record_set):
        if filterset['nonuniq'][1] == 'highest-GenTrain':
            maxval = max(rec.INFO.get('GenTrain_Score') for rec in record_set)
            attr_get = lambda rec: rec.INFO.get('GenTrain_Score')
        elif filterset['nonuniq'][1] == 'highest-GenCall':
            maxval = max(c.data.get('IGC', 0) for rec in record_set for c in rec.calls)
            attr_get = lambda rec: max(c.data.get('IGC', 0) for c in rec.calls)
        else:
            maxval = None
            attr_get = lambda rec: None

        for rec in record_set:
            if attr_get(rec)== maxval:
                continue
            rec.add_filter(filterset['nonuniq'][0])

    prev_pos = (None, None)
    record_set = []
    for record in reader:
        # Apply GT &n (I)GC filters
        # import pdb; pdb.set_trace()
        if 'GT' in filterset:
            if record.INFO.get('GenTrain_Score') < filterset['GT'][1]:
                record.add_filter(filterset['GT'][0])
        if 'GC' in filterset:
            gc_scores = (c.data.get('IGC', 0) for c in record.calls)
            if any(score < filterset['GC'][1] for score in gc_scores):
                record.add_filter(filterset['GC'][0])

        # Assume records are sorted by position
        # If same position, wait (recored is saved in record_set)
        if (record.POS, record.CHROM) == prev_pos:
            logger.debug(f"Skipping write on {record.CHROM}:{record.POS}, {record.ID}")
        # If new position (except first), apply uniq_pos filter & write out _previous_ records
        elif prev_pos != (None, None):
            if 'nonuniq' in filterset:
                apply_uniq_pos_filter(record_set)
            logger.debug(f"Writing for {record_set[0].CHROM}:{record_set[0].POS} with {len(record_set)} records")
            for rec in record_set:
                writer.write_record(rec)
            record_set = []
        # Add current record to set (same pos or just emptied)
        record_set.append(record)
        prev_pos = (record.POS, record.CHROM)

    # Write out last record set
    if 'nonuniq' in filterset:
        apply_uniq_pos_filter(record_set)
    for rec in record_set:
        writer.write_record(rec)

    writer.close()



filterset = get_filter_settings(snakemake.wildcards['filter'], snakemake.config)

for name, vals in filterset.items():
    logger.info(f"Adding Filter for {name} as '{vals[0]}'")

filter_snp_vcf(filterset, snakemake.input[0], snakemake.output[0])
