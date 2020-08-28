## MySQL configuration notes for gEAR instances

This document is intended to capture all the customizations to the mysql config files needed for a gEAR instance to operate properly.  

### Installation

    $ sudo apt install mysql-server

### Setting a root password

This isn't entirely required, but I had to chase this issue down for some time.  On Ubuntu 18.04 and higher the root user is managed via the auth_socket plugin instead of mysql_native_password.  To fix this:


    mysql> ALTER USER 'root'@'localhost' IDENTIFIED WITH mysql_native_password BY 'whatever';
    Query OK, 0 rows affected (0.00 sec)

    mysql> flush privileges;
    Query OK, 0 rows affected (0.00 sec)


