$ sudo mount -o discard,defaults /dev/sdb /mnt/disks/epiviz-drive

### RabbitMQ consumer startup
sudo su -

cd ~jorvis/git/gEAR/listeners/

nohup /opt/bin/python3 projectr_consumer.py >>/var/log/gEAR_queue/projectr.log 1>/dev/null 2>>/var/log/gEAR_queue/projectr.log &
