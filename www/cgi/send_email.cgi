#!/opt/bin/python3

"""
Sends email to user
Scopes: 'forgot_password'

"""

import cgi, json
import configparser
import hashlib
import sys
import mysql.connector
from mysql.connector import errorcode
import uuid

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

def main():

    print('Content-Type: application/json\n\n')
    result = {'error': [], 'success': 0 }

    config = configparser.ConfigParser()
    config.read('../../gear.ini')

    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'])
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password", file=sys.stderr)
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist", file=sys.stderr)
        else:
            print(err, file=sys.stderr)

    cursor = cnx.cursor()
    form = cgi.FieldStorage()
    email = form.getvalue('email')
    scope = form.getvalue('scope')
    destination_page = form.getvalue('destination_page')

    if scope == 'forgot_password':

        user_help_id = get_help_id(cursor, email) # return 'False' or help_id string
        if user_help_id:
            #email was found.

            # url to change password
            url = destination_page + '?help_id=' + user_help_id

            # https://docs.python.org/3/library/email-examples.html
            gear = 'gearportal.igs@gmail.com'
            # user = email
            msg = MIMEMultipart('alternative')
            msg['Subject'] = 'Reset gEAR Password'
            msg['From'] = gear
            msg['To'] = email

            text = "You requested to change your gEAR password. Please click the link to continue:  {0}".format(url)
            html = """\
            <html>
            <head></head>
            <body style="font-family:Helvetica Neue,Helvetica,Arial,sans-serif;">
            <div style="height:50px; background-color:#2F103E;">
               <a href="http://umgear.org/index.html">
                 <img src="http://umgear.org/img/logo_h70_w126.png" style="border-radius:4px; margin-left:90px;">
              </a>
            </div>
            <div style="text-align:center; height:50vh; vertical-align:middle; margin-top:100px;">
              <br />
              <h2 style="font-weight:bold;color: rgb(134,134,134);">Forgot your password?</h2>
              <p> Click below to create a new one.</p>
              <a href='{0}'>Create new password</a>
              <br /><br />
              <p>Please do not reply to this message.</p>
            </div>
            </body>
            </html>
            """.format(url)

            part1 = MIMEText(text, 'plain')
            part2 = MIMEText(html, 'html')

            msg.attach(part1)
            msg.attach(part2)

            # http://stackoverflow.com/a/17596848/2900840
            s = smtplib.SMTP('smtp.gmail.com:587')
            s.ehlo()
            s.starttls()
            s.login('gearportal.igs@gmail.com', 'gearadmin')

            s.sendmail( gear, email, msg.as_string() )
            s.quit()

            result['success'] = 1
        else:
            #email was not in database.
            result['error'] = "We could not find that email. Please check what you entered and try again."

    print(json.dumps(result))

    cnx.commit()
    cursor.close()
    cnx.close()


def get_help_id(cursor, email):
    help_id = None

    qry = "SELECT help_id FROM guser WHERE email = %s"
    cursor.execute(qry, (email,))
    for row in cursor:
        help_id = row[0]

    return help_id


if __name__ == '__main__':
    main()
