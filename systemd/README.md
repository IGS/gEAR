# system.d files

These are service files to place in `/etc/systemd/system/` to run as a service.

The "@" in the service file name indicates that this file is a template file.  You can spawn off workers off of this template by appending numbers after the "@". For instance, if you have a "worker@.service" template, you can spawn off "worker@2.service", "worker@2.service", etc. You can also just use the ".target" file to group all these workers together as a service rather than having to start each indiviudally. See https://www.stevenrombauts.be/2019/01/run-multiple-instances-of-the-same-systemd-unit/ for more info on all this.

## Non-persisting start

To start the service run `sudo systemctl start <service_filename>`. If any changes are made to the service run `sudo systemctl daemon-reload`, and then run the "start" command again. This will not persist during a reboot.

## Reboot-persisting start

To start a service that persists during a reboot, run `sudo systemctl enable <service_filename>`

## Notes

* For a service that is supposed to persist, like a RabbitMQ consumer, use Service->Type=Simple.  If you use Type->Forking, I believe systemctl will forever wait for the forked process to exit and return to command line.