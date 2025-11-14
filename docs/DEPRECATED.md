# Deprecated Documentation Files

This document lists files in this directory that have been superseded by the new developer documentation structure.

## Migration to `/docs/developer/`

The developer documentation has been reorganized into a more coherent structure at `/docs/developer/`. Please use the new documentation instead of the files listed below.

## Setup Guides (Moved)

The following setup guides have been **moved** to `/docs/developer/setup/`:

- ~~`setup.new_server.notes.md`~~ → [`docs/developer/setup/new_server.md`](developer/setup/new_server.md)
- ~~`setup.mysql.md`~~ → [`docs/developer/setup/mysql.md`](developer/setup/mysql.md)
- ~~`setup.apache.md`~~ → [`docs/developer/setup/apache.md`](developer/setup/apache.md)
- ~~`setup.python.md`~~ → [`docs/developer/setup/python.md`](developer/setup/python.md)
- ~~`setup.r_rpy2.md`~~ → [`docs/developer/setup/r_rpy2.md`](developer/setup/r_rpy2.md)
- ~~`setup.rabbitmq.md`~~ → [`docs/developer/setup/rabbitmq.md`](developer/setup/rabbitmq.md)
- ~~`setup.epiviz.md`~~ → [`docs/developer/setup/epiviz.md`](developer/setup/epiviz.md) (**OBSOLETE** - Epiviz removed, replaced by Gosling)

See the [Setup Guide Index](developer/setup/README.md) for the complete setup documentation.

## Component Documentation (Moved)

The following component-specific docs have been integrated into the new structure:

### Docker
- Original: `docker/docker_notes.md`
- New location: [`docs/developer/setup/docker.md`](developer/setup/docker.md)

### Systemd Services
- Original: `systemd/README.md`
- New location: [`docs/developer/setup/systemd.md`](developer/setup/systemd.md)

### Microservices
- Original: `services/projectr/README.md`
- New location: [`docs/developer/services/projectr.md`](developer/services/projectr.md)
- See also: [`docs/developer/services/README.md`](developer/services/README.md) for all services

## Legacy Files (Likely Obsolete)

The following files may be outdated and are candidates for removal pending team review:

### Obsolete/Removed
- `epiviz_notes-old.md` - Old Epiviz notes (**OBSOLETE** - Epiviz has been removed and replaced by Gosling)
- `epiviz-tracks.json` - Epiviz configuration (**OBSOLETE** - Epiviz removed)
- `projectR_scratch_notes.md` - Scratch notes (superseded by formal ProjectR docs)

### Specialized/Internal Documentation
These files contain specialized knowledge and should be reviewed before removal:

- `adding_new_display_types.md` - Developer guide for adding new display types
- `copying_gene_symbols.md` - Internal process documentation
- `cron_config.md` - Cron job configuration
- `gene_curator_notes.md` - Curator-specific notes
- `how_analyses_panels_are_displayed.md` - Technical documentation
- `multigene_curator.md` - Multi-gene curator documentation
- `mysql_config.md` - MySQL configuration notes
- `plugins.md` - Plugin system documentation
- `release_test_plan.md` - Release testing procedures
- `svg_formatting.md` - SVG formatting guidelines
- `uploading_spatial_dataset.md` - Spatial dataset upload process

**Action Required**: These specialized docs should be:
1. Reviewed for current relevance
2. Integrated into the new developer documentation where appropriate
3. Moved to a "technical-notes" subdirectory if still relevant but not fitting the new structure
4. Removed if obsolete

## WWW Docs Directory

The `/www/docs/` directory contains **end-user documentation** for the gEAR wiki (Docsify-based).

**Important**: Do NOT modify files in this directory unless updating end-user documentation.

## Recommendations

### For Developers
1. Use [`docs/developer/README.md`](developer/README.md) as your starting point
2. Refer to the new structured documentation in `docs/developer/`
3. If you find information missing, please update the new docs or create an issue

### For Maintainers (@adkinsrs)
1. Review specialized documentation files listed above
2. Determine which should be:
   - Integrated into new developer docs
   - Kept as-is in a technical notes directory
   - Removed as obsolete
3. Consider creating a `docs/technical-notes/` directory for internal process docs

### File Removal Timeline
- **Phase 1** (Immediate): Add deprecation notices to moved setup files
- **Phase 2** (After team review): Move/integrate specialized docs
- **Phase 3** (After validation): Remove obsolete files

## Questions?

If you're unsure whether a file is still needed:
1. Check the new developer documentation first
2. Ask in team channels or check with @adkinsrs
3. Create an issue on GitHub for clarification
