## Cron config for gEAR instances

This document is intended to list all requires for setting up cron for gEAR instances.

### Setting up the cron entries

    $ crontab -e

    */30 * * * * cd /home/jorvis/git/gEAR/bin && ./generate_initial_composition_plots.py

### Cron logs.  Can depend on the system, but you can find them like this:

    $ grep -ic cron /var/log/* | grep -v :0
     ./cron:52
     ./dracut.log:2
     ./messages:1
     ./secure:6



