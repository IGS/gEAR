# gEAR Documentation Index

Welcome to the gEAR documentation. This index will help you find the right documentation for your needs.

## For End Users

**gEAR User Documentation** → [GitHub Wiki](https://github.com/IGS/gEAR/wiki) or visit any gEAR portal's [manual page](https://umgear.org/manual.html)

Topics covered:
- Getting started with gEAR
- Searching for genes and datasets
- Creating and managing dataset collections
- Using analysis tools
- Uploading and curating data
- Video tutorials and presentations

**Do not modify files in `docs/Documentation/`** - These are synced to the wiki.

## For Developers

**Main Developer Guide** → [`developer/README.md`](developer/README.md)

Comprehensive guide covering:
- Architecture overview
- Development workflow
- Code style guidelines
- Testing procedures
- Common development tasks

### Quick Links for Developers

- **[Setup Guides](developer/setup/README.md)** - Server installation and configuration
  - New server setup
  - MySQL, Apache, Python, R, RabbitMQ
  - Docker development environment
  - Systemd services

- **[Utility Scripts](developer/scripts/README.md)** - Documentation for 105+ scripts in `/bin`
  - Data conversion and format scripts
  - H5AD manipulation tools
  - Database loading utilities
  - Validation and testing scripts
  - Dataset management tools

- **[Microservices](developer/services/README.md)** - Service documentation
  - ProjectR service (matrix projection)
  - Spatial panel service (spatial transcriptomics)
  - RabbitMQ consumers (background workers)

## Documentation Organization

```
docs/
├── README.md (this file)          # Documentation index
├── developer/                     # Developer documentation
│   ├── README.md                  # Main developer guide
│   ├── setup/                     # Server setup guides
│   ├── scripts/                   # Bin scripts documentation
│   └── services/                  # Microservices documentation
├── Documentation/                 # End-user wiki (DO NOT MODIFY)
├── posters/                       # Historical presentations (DO NOT MODIFY)
├── ui-v2-design/                  # UI v2 prototypes (DO NOT MODIFY)
├── DEPRECATED.md                  # List of superseded files
└── [various .md files]            # Legacy/specialized docs
```

## Legacy Documentation

Some files in the main `docs/` directory have been superseded by the new developer documentation structure. See [`DEPRECATED.md`](DEPRECATED.md) for details on which files have been moved or are obsolete.


**Note**: These files are under review for integration into the new structure. Use with caution as some may be outdated.

## Component-Specific Documentation

Some directories have their own documentation files with notes specific to that component:

- `docker/docker_notes.md` - Docker setup notes (now at `developer/setup/docker.md`)
- `systemd/README.md` - Systemd service information (now at `developer/setup/systemd.md`)
- `services/projectr/README.md` - ProjectR service (now at `developer/services/projectr.md`)

These original files are kept for reference but now include notices pointing to the consolidated documentation.

## Finding What You Need

### "I want to set up my own gEAR instance"
→ Start with [`developer/setup/README.md`](developer/setup/README.md)

### "I want to understand how gEAR works"
→ Read [`developer/README.md`](developer/README.md) - Architecture Overview section

### "I need to run a specific script from /bin"
→ Check [`developer/scripts/README.md`](developer/scripts/README.md) for documentation

### "I'm having issues with ProjectR/Spatial/RabbitMQ services"
→ See [`developer/services/README.md`](developer/services/README.md)

### "I want to contribute to gEAR"
→ Start with [`developer/README.md`](developer/README.md) - Development Workflow section

### "I'm a gEAR user looking for help"
→ Visit the [GitHub Wiki](https://github.com/IGS/gEAR/wiki) or [UMgEAR manual](https://umgear.org/manual.html)

## Contributing to Documentation

If you find missing, outdated, or incorrect documentation:

1. **For developer docs**: Update files in `docs/developer/` and submit a PR
2. **For end-user docs**: Update the [GitHub Wiki](https://github.com/IGS/gEAR/wiki)
3. **For questions**: Create an issue on GitHub or ask in team channels

### Documentation Standards

When adding or updating documentation:

- Use clear, descriptive headings
- Include code examples where appropriate
- Link to related documentation
- Keep information current (note last update date for time-sensitive content)
- Follow existing formatting and style

## Additional Resources

- **GitHub Repository**: https://github.com/IGS/gEAR
- **Issue Tracker**: https://github.com/IGS/gEAR/issues
- **Production Instances**:
  - [UMgEAR](https://umgear.org) - Hearing research portal
  - [NeMO Analytics](https://nemoanalytics.org) - Brain research portal

## Getting Help

1. Check the appropriate documentation section above
2. Search the [GitHub Wiki](https://github.com/IGS/gEAR/wiki) for end-user questions
3. Review existing [GitHub Issues](https://github.com/IGS/gEAR/issues)
4. Ask in team channels (for team members)
5. Create a new GitHub issue with appropriate labels

---

**Last Updated**: October 2025
**Maintainer**: @adkinsrs maintains infrastructure documentation
