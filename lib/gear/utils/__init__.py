"""
gear.utils - Cross-cutting helper modules shared across multiple scripts, grouped by concern:

- gene_mapping: AnnData/var gene-symbol -> Ensembl ID mapping
- obs: obs-dataframe sanitization/categorization helpers
- resource_limits: process memory-limit management
- job_coordination: RabbitMQ consumer job-lock/retry/logging helpers
- fulltext: user search text -> MySQL full-text boolean-mode query
- archives: user-facing errors for uploaded archives that can't be extracted

Import from the specific submodule you need, e.g.:
    from gear.utils.job_coordination import log_line
"""
