from datetime import datetime
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage

gear_email = 'gearportal.igs@gmail.com'
gear_email_pass = 'gearadmin'

# base template for all gEAR emails
base_html = """\
    <html>
    <body style="font-family:Helvetica Neue,Helvetica,Arial,sans-serif;">
    <div style="height:50px; background-color:#2F103E;">
       <a href="https://umgear.org/index.html">
         <img src="https://umgear.org/img/logo_h70_w126.png" style="border-radius:4px; margin-left:90px;">
      </a>
    </div>
    <div style="text-align:center; height:50vh; vertical-align:middle; margin-top:100px;">
      {html_message}
      <br /><br />
      <p>Please do not reply to this message.</p>
    </div>
    </body>
    </html>
"""

# Message template: load dataset
html_message_upload = """\
    <h2 style="font-weight:bold;color:rgb(134,134,134);">Dataset Upload {load_status_titled}</h2>
    <br />
    <p>Your upload for dataset: <b>{dataset_title}</b> has <b>{load_status}</b>.</p>
    <br />
    <p style="color:red;"><b>{error_message}</b></p>
    <br />
    <p>Click below to go to the Dataset Manager:</p>
    <a href='https://umgear.org/dataset_manager.html' target='_blank'>Dataset Manager</a>
    <br /><br />
    <p>For help, please visit the gEAR <a href='https://umgear.org/manual.html' target='_blank'> gEAR manual</a> or <a href='https://umgear.org/contact.html' target='_blank'>contact us</a>.</p>
"""

class Email:
    def __init__(self, sender_email=None, sender_password=None, receiver_email=None, subject=None, text_message=None, html_message=None):
        self.sender_email = sender_email
        self.sender_password = sender_password
        self.receiver_email = receiver_email
        self.subject = subject
        self.text_message = text_message
        self.html_message = html_message

        #Default is gEAR email and password
        if self.sender_email is None:
            self.sender_email = gear_email
        if self.sender_password is None:
            self.sender_password = gear_email_pass

    def send(self):
        '''
        Sends the email to the recipient.

        First builds the email message, then connects to gmail and sends email.

        Note: Only Gmail email addresses will work. Also, if using a sender_email
            other than the gEAR address may error due to that account's security settings (i.e 2-factor authentication)
        '''

        #Build email message
        message = self._build_message()

        # http://stackoverflow.com/a/17596848/2900840
        s = smtplib.SMTP('smtp.gmail.com:587')
        s.ehlo()
        s.starttls()
        try:
            s.login(self.sender_email, self.sender_password)
            s.sendmail( self.sender_email, self.receiver_email, message.as_string() )
            s.quit()
        except Exception as err:
            error = str(err)
            raise Exception("{0}\t{1}: Unable to send email: {2}".format(error))


    def _build_message(self):
        '''
        Builds message component of email.
        Requires:
            sender_email, receiver_email, subject, text_message, html_message
        '''
        if self.sender_email is None:
            raise Exception("Cannot build email message. No 'sender_email' found")
        if self.receiver_email is None:
            raise Exception("Cannot build email message. No 'receiver_email' found")
        if self.subject is None:
            raise Exception("Cannot build email message. No 'subject' found")
        if self.text_message is None:
            raise Exception("Cannot build email message. No 'text_message' found")
        if self.html_message is None:
            raise Exception("Cannot build email message. No 'html_message' found")

        # https://docs.python.org/3/library/email-examples.html
        msg = MIMEMultipart('alternative')
        msg['From'] = self.sender_email
        msg['To'] = self.receiver_email
        msg['Subject'] = self.subject
        msg.attach(MIMEText(self.text_message, 'plain'))
        msg.attach(MIMEText(self.html_message, 'html'))

        # Add gEAR logo
        fp = open('../www/img/logo_h70_w126.png', 'rb')
        msg_image = MIMEImage(fp.read())
        fp.close()
        msg_image.add_header('Content-ID', '<gear_logo>') #<img src="cid:gear_logo">
        msg.attach(msg_image)

        return msg
