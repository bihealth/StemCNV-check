import vcfpy
from loguru import logger
from collections import defaultdict

logger.remove()
# use at least INFO level here
logger.add(snakemake.log['err'], level=["WARNING", "INFO", "DEBUG"][min(snakemake.config['verbose_level']+1, 2)])
logger.info(f"Reading in '{snakemake.input[0]}' and will write to '{snakemake.output[0]}'")


def get_filter_settings(filtersetname, config, sample_sex):

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

    pseudoauto = config_filterset.get('Pseudoautosomal')
    if pseudoauto == 'remove' or (pseudoauto == 'remove-male' and sample_sex == 'm'):
        filterset['pseudoauto'] = ('PSEUDO-AUTOSOMAL', pseudoauto, 'Remove records in pseudoautosomal regions')

    if sample_sex == 'f':
        filterset['FEMALE-Y'] = ('FEMALE-Y', None, 'Remove Y records in female sample')

    return filterset


def apply_uniq_pos_filter(record_set):
    logger.debug(f"Applying non-unique position filter to {len(record_set)} records at {record_set[0].CHROM}:{record_set[0].POS}")
    if FILTERSET['nonuniq'][1] == 'highest-GenTrain':
        attr_get = lambda rec: rec.INFO.get('GenTrain_Score')
    elif FILTERSET['nonuniq'][1] == 'highest-GenCall':
        attr_get = lambda rec: rec.calls[0].data.get('IGC', 0)
    else:
        return
    # records themselves are comparable, disfavour './.' entries
    # (which are either split multi-allelic or otherwise less informative)
    rec_selected = max(
        record_set,
        key=lambda rec: (
            attr_get(rec) -
            1 if rec.calls[0].data.get('GT', 0) == './.' else 0
        )
    )
    # Add filter to all but the selected record
    for rec in record_set:
        if rec != rec_selected:
            rec.add_filter(FILTERSET['nonuniq'][0])
        else:
            logger.debug(f"Not adding filter to selected record {rec.ID[0]} with '{FILTERSET['nonuniq'][1]}' = {attr_get(rec)}")


def apply_pseudo_autosomal_filter(record):
    if record.CHROM[-1] in ('X', 'Y'):
        for start, end in PSEUDOAUTO_REGIONS[GENOME_VERSION][record.CHROM[-1]]:
            if start <= record.POS <= end:
                record.add_filter(FILTERSET['pseudoauto'][0])
                break





def filter_snp_vcf(filterset, sample_sex, input='/dev/stdin', output='/dev/stdout'):

    # Function to handle/finalize set of recors with the same position
    def finalize_record_set(record_set):
        if 'nonuniq' in filterset and len(record_set) > 1:
            apply_uniq_pos_filter(record_set)

        logger.debug(f"Writing for {record_set[0].CHROM}:{record_set[0].POS} with {len(record_set)} records")
        # ID's need to be unique, but splitting multi-allelic sites introduces duplicates
        # -> need to modify ID's to be unique
        rec_ids = defaultdict(int)
        for rec in record_set:
            rec_ids[rec.ID[0]] += 1
            if rec_ids[rec.ID[0]] > 1:
                rec.ID[0] = rec.ID[0] + f'_{rec_ids[rec.ID[0]]}'
            if not rec.FILTER:
                rec.add_filter('PASS')
            writer.write_record(rec)

    with open(input) as input_stream:
        reader = vcfpy.Reader.from_stream(input_stream)

        for _, filtercond in filterset.items():
            reader.header.add_filter_line(
                vcfpy.OrderedDict([
                    ('ID', filtercond[0]),
                    ('Description', filtercond[2])
                ])
            )

        with open(output, 'w') as output_stream:
            writer = vcfpy.Writer.from_stream(output_stream, reader.header)
            prev_pos = (None, None)
            record_set = []
            for record in reader:
                # Fail on multi-sample vcf
                if len(record.calls) > 1:
                    raise ValueError("Multi-sample VCF not supported")
                # Always hard-filter records without proper REF definition
                # these are mainly Indel-Probes with "I" or "D" in REF & ALT
                if record.REF not in ['A', 'C', 'G', 'T']:
                    continue

                # Always filter female Y probes
                if sample_sex == 'f' and record.CHROM[-1] == 'Y':
                    record.add_filter(FILTERSET['FEMALE-Y'][0])

                # apply pseudo-autosomal if turned on (=in filterset) and:
                # used for all samples ('remove') or (implicit last option) only matching male sample
                if 'pseudoauto' in filterset and (filterset['pseudoauto'][1] == 'remove' or sample_sex == 'm'):
                    apply_pseudo_autosomal_filter(record)

                # Apply GT & (I)GC filters
                if 'GT' in filterset and record.INFO.get('GenTrain_Score') < filterset['GT'][1]:
                    record.add_filter(filterset['GT'][0])
                if 'GC' in filterset and record.calls[0].data.get('IGC', 0) < filterset['GC'][1]:
                    record.add_filter(filterset['GC'][0])
                # Assume records are sorted by position
                # If same position, wait (record is saved in record_set)
                if (record.POS, record.CHROM) == prev_pos:
                    logger.debug(f"Skipping write on {record.CHROM}:{record.POS}, {record.ID[0]}")
                # If new position (except first), apply uniq_pos filter & write out _previous_ records
                elif prev_pos != (None, None):
                    finalize_record_set(record_set)
                    record_set = []

                # Add current record to set (same pos or just emptied)
                record_set.append(record)
                prev_pos = (record.POS, record.CHROM)

            # Write out last record set
            finalize_record_set(record_set)


PSEUDOAUTO_REGIONS = {
    "hg38": {
        'X': [(10001, 2781479), (155701383, 156030895)],
        'Y': [(10001, 2781479), (56887903, 57217415)],
    },
    "hg19": {
        'X': [(60001, 2699520), (88484850, 92327352), (154931044, 155260560)],
        'Y': [(10001, 2649520), (59034050, 59363566)],
    }
}

FILTERSET = get_filter_settings(snakemake.wildcards['filter'], snakemake.config, snakemake.params['sample_sex'])
GENOME_VERSION = 'hg38' if snakemake.params['genome_version'] in ('hg38', 'GRCh38') else 'hg19'

for name, vals in FILTERSET.items():
    logger.info(f"Adding Filter for {name} as '{vals[0]}'")

filter_snp_vcf(FILTERSET, snakemake.params['sample_sex'], snakemake.input[0], snakemake.output[0])
