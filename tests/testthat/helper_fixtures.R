suppressMessages(library(tidyverse))

minimal_probes <- tibble(
    list_name = 'test-list',
    hotspot = c('DDX11L1', 'dummyC', 'chr1:40000-50000', '1p36', '1p35.2', '1q21'),
    mapping = c('gene_name', 'gene_name', 'position', 'gband', 'gband', 'gband'),
    call_type = c('any', 'gain', 'any', 'loss', 'any', 'gain'),
    check_score = c(30, 15, 10, 10, 12, 10),
    description = c('Sources: dummy', 'Something: Dummy{1}', 'dummy', 'Sources: dummy{1}\\nSomething: else{2}',
                    'dummy', 'Sources: dummy{1},dummy{2}'),
    description_doi = c(NA, 'link', NA, 'link, link2', NA, 'link, link2'),
    description_htmllinks = c(
        'Sources: dummy', 'Something: <a href="link">Dummy</a>', 'dummy',
        'Sources: <a href="link">dummy</a>\\nSomething: <a href="link2">else</a>',
        'dummy', 'Sources: <a href="link">dummy</a>,<a href="link2">dummy</a>'
    )
)