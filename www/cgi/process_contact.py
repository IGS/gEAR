#!/opt/bin/python3

import cgi
import datetime

form = cgi.FieldStorage()
submitter_name = form.getfirst('submitterName')

name = "{0}".format(datetime.datetime.now())
name = name.replace(' ', '_')

ofh = open("./user-submitted-comments/{0}".format(name), 'wt')

ofh.write("submitter: {0}\n".format(submitter_name))
ofh.write("e-mail: {0}\n".format(form.getfirst('InputEmail')))
ofh.write("keep updated: {0}\n".format(form.getfirst('keep_updated')))
ofh.write("sec check: {0}\n".format(form.getfirst('super_impressive_security_check')))
ofh.write("comments: {0}\n".format(form.getfirst('comments')))

print("Content-type: text/html\n")
print("Comment submitted.<br><br>Back to <a href='../index.html'>main site</a>.")

