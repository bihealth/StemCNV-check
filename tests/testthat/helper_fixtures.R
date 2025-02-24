# library(tidyverse)

minimal_probes <- tibble::tibble(
    list_name = 'test-list',
    hotspot = c('DDX11L1', 'dummyC', 'chr1:40000-50000', '1p36', '1p35.2', '1q21'),
    mapping = c('gene_name', 'gene_name', 'position', 'gband', 'gband', 'gband'),
    call_type = c('any', 'gain', 'any', 'loss', 'any', 'gain'),
    check_score = c(30, 15, 10, 10, 12, 10),
    description = c('Sources: dummy', 'Something: Dummy{1}', 'dummy', 'Sources: dummy{1}\\nSomething: else{2}',
                    'dummy', 'Sources: dummy{1},dummy{2}'),
    description_doi = c(NA, 'doi', NA, 'doi, doi2', NA, 'doi, doi2'),
    description_htmllinks = c(
        'Sources: dummy', 
        'Something: <a href="https://doi.org/doi" target="_blank" rel="noopener noreferrer">Dummy</a>',
        'dummy',
        'Sources: <a href="https://doi.org/doi" target="_blank" rel="noopener noreferrer">dummy</a>&#013;Something: <a href="https://doi.org/doi2" target="_blank" rel="noopener noreferrer">else</a>',
        'dummy', 
        'Sources: <a href="https://doi.org/doi" target="_blank" rel="noopener noreferrer">dummy</a>,<a href="https://doi.org/doi2" target="_blank" rel="noopener noreferrer">dummy</a>'
    )
)