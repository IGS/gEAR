$ sudo mount -o discard,defaults /dev/sdb /mnt/disks/epiviz-drive

### RabbitMQ consumer startup
sudo su -

## SAdkins - commenting out since service files are available
#cd ~jorvis/git/gEAR/listeners/
#nohup /opt/bin/python3 projectr_consumer.py >>/var/log/gEAR_queue/projectr.log 1>/dev/null 2>>/var/log/gEAR_queue/projectr.log &

Go to /etc/systemd/system/projectr-consumer.service. Make sure the projectr daemon is set to "restart=always". If not, fix it and run `sudo systemctl daemon-reload`

sudo systemctl enable projectr-consumer
sudo systemctl start projectr-consumer