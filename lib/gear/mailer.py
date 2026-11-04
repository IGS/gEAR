"""
mailer.py - Sending e-mail from gEAR.

Uses the account in gear.ini [email_sender] (address, password, and optionally smtp_server,
default smtp.gmail.com:587). Curator notifications go to the comma-separated addresses in
gear.ini [curation] notify_emails.
"""

import smtplib
import sys
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from gear.serverconfig import ServerConfig

DEFAULT_SMTP_SERVER = "smtp.gmail.com:587"


def send_email(recipients: list[str], subject: str, text: str, html: str | None = None) -> None:
    """
    Send an e-mail with a plain text body and, optionally, an HTML alternative.

    Raises an exception if the server config is missing or delivery fails.
    """
    config = ServerConfig().parse()
    sender = config["email_sender"]["address"]
    password = config["email_sender"]["password"]
    smtp_server = config.get("email_sender", "smtp_server", fallback=DEFAULT_SMTP_SERVER)

    # https://docs.python.org/3/library/email-examples.html
    msg = MIMEMultipart("alternative")
    msg["From"] = sender
    msg["To"] = ", ".join(recipients)
    msg["Subject"] = subject
    msg.attach(MIMEText(text, "plain"))
    if html is not None:
        msg.attach(MIMEText(html, "html"))

    # http://stackoverflow.com/a/17596848/2900840
    s = smtplib.SMTP(smtp_server)
    s.ehlo()
    s.starttls()
    s.login(sender, password)
    s.sendmail(sender, recipients, msg.as_string())
    s.quit()


def get_curator_emails() -> list[str]:
    """Addresses in gear.ini [curation] notify_emails."""
    config = ServerConfig().parse()
    value = config.get("curation", "notify_emails", fallback="")
    return [address.strip() for address in value.split(",") if address.strip()]


def notify_curators(subject: str, text: str) -> bool:
    """
    E-mail the gEAR curators. Returns False, after logging why, if no curators are configured
    or delivery fails, so callers can carry on.
    """
    recipients = get_curator_emails()
    if not recipients:
        print(f"No [curation] notify_emails configured; not sending: {subject}", file=sys.stderr)
        return False

    try:
        send_email(recipients, subject, text)
    except Exception as e:
        print(f"Could not e-mail curators ({subject}): {e}", file=sys.stderr)
        return False
    return True
