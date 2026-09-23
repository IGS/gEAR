# Moved and Removed Documentation

Documentation that used to live directly in `docs/` or in component directories has been moved or removed. Start from the [documentation index](README.md) or the [developer guide](developer/README.md).

## Moved

| Old location | New location |
| --- | --- |
| `docs/setup.new_server.notes.md` | [developer/setup/new_server.md](developer/setup/new_server.md) |
| `docs/setup.mysql.md` | [developer/setup/mysql.md](developer/setup/mysql.md) |
| `docs/setup.apache.md` | [developer/setup/apache.md](developer/setup/apache.md) |
| `docs/setup.python.md` | [developer/setup/python.md](developer/setup/python.md) |
| `docs/setup.r_rpy2.md` | [developer/setup/r_rpy2.md](developer/setup/r_rpy2.md) |
| `docs/setup.rabbitmq.md` | [developer/setup/rabbitmq.md](developer/setup/rabbitmq.md) and [developer/services/rabbitmq_consumers.md](developer/services/rabbitmq_consumers.md) |
| `docker/docker_notes.md` | [developer/setup/docker.md](developer/setup/docker.md) |
| `docker/mysql_setup_notes.md` | [developer/setup/docker_mysql.md](developer/setup/docker_mysql.md) |
| `systemd/README.md` | [developer/setup/systemd.md](developer/setup/systemd.md) |
| `services/projectr/README.md` | [developer/services/projectr.md](developer/services/projectr.md) |
| `docs/plugins.md` | [developer/services/plugins.md](developer/services/plugins.md) |
| `docs/webpage_dependencies.md` | [developer/misc/webpage_dependencies.md](developer/misc/webpage_dependencies.md) |
| `docs/release_test_plan.md` | [developer/misc/release_test_plan.md](developer/misc/release_test_plan.md) |
| `docs/adding_new_display_types.md` | [misc/adding_new_display_types.md](misc/adding_new_display_types.md) |
| `docs/copying_gene_symbols.md` | [misc/copying_gene_symbols.md](misc/copying_gene_symbols.md) |
| `docs/gene_curator_notes.md` | [analyst/gene_curator_notes.md](analyst/gene_curator_notes.md) |
| `docs/how_analyses_panels_are_displayed.md` | [analyst/how_analyses_panels_are_displayed.md](analyst/how_analyses_panels_are_displayed.md) (historical) |
| `docs/multigene_curator.md` | [analyst/multigene_curator.md](analyst/multigene_curator.md) |
| `docs/svg_formatting.md` | [analyst/svg_formatting.md](analyst/svg_formatting.md) |
| `docs/uploading_spatial_dataset.md` | [analyst/uploading_spatial_dataset.md](analyst/uploading_spatial_dataset.md) |

## Removed

- `docs/mysql_config.md` - MySQL configuration notes, merged into [developer/setup/mysql.md](developer/setup/mysql.md)
- `docs/cron_config.md` - cron setup notes, removed
- `docs/projectR_scratch_notes.md` - superseded by [developer/services/projectr.md](developer/services/projectr.md)
- `docs/epiviz_notes-old.md`, `docs/epiviz-tracks.json`, `docs/developer/setup/epiviz.md` - Epiviz was removed and replaced by Gosling
- `docs/REORGANIZATION_SUMMARY.md` - one-off summary of the documentation reorganization
- `www/manual.html` (with `www/js/manual.js`, `www/css/manual.css` and `www/include/manual/`) - legacy v1 user manual page, replaced by the [wiki](wiki/gEARWiki.md) and the video guides on the home page
