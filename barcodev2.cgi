#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#Selection of barcodes for GBS v2.0 written by Thomas p van Gurp
import cgitb
import subprocess
cgitb.enable()


import cgi
import re
import sys
import sys,random,Bio
sys.path.append("/modules")
sys.path.append(".")
import smtplib
from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from optparse import OptionParser
import csv
import time, datetime
print "time.localtime() : %s" % time.localtime()

form = cgi.FieldStorage()

def email_success(opts):
    """Send an email that barcodes are succesfully generated"""
    output_location = "<a href='%s'>Download markers</a>"%(opts.output.replace('/httpdocs',''))
    #print 'Download csv-output here: %s'%(output_location)
    sender = 'info@deenabio.com'
    receivers = ['%s'%opts.email, 'thomas@deenabio.com']
    #print opts
    message = """From: Deena Bioinformatics <info@deenabio.com>\nTo:  %s <%s>\nSubject: Your %s-plex %s GBS Barcodes
    \nMIME-Version: 1.0
    \nContent-type: text/html

    Dear %s,
    <br><br>
    Thanks for using my Barcode generating script! 
    <br><br>
    You can now download your Barcodes at: <a href='http:\\www.deenabio.com/%s'>www.deenabio.com/%s</a>
    <br><br>
    with kind regards,
    <br><br>
    <a href='http://nl.linkedin.com/in/thomasvangurp'>Thomas van Gurp</a><br>
    Bioinformatics consultant | Owner<br>
    Deena Bioinformatics<br>
    <a href='http:\\www.deenabio.com'>www.deenabio.com</a>
    """%(opts.name, opts.email, opts.number, opts.enzyme, opts.name,\
         opts.output.replace('/httpdocs','')[3:],opts.output.replace('/httpdocs','')[3:])
    try:
       smtpObj = smtplib.SMTP('localhost')
       smtpObj.sendmail(sender, receivers[0], message)         
       print """A download link for the Barcodes was sent to you by email
       
       Thanks for your patience!"""
    except SMTPException:
       print """Error: unable to send email
       
       Please provide a valid email-adress"""
    smtpObj.sendmail(sender, receivers[1], message)

def email_fail(opts, error):
    """Send an error report that barcodes were not succesfully generated"""
    sender = 'info@deenabio.com'
    receivers = ['thomas@deenabio.com']
    message = """From: Deena Bioinformatics <info@deenabio.com>\nTo:  %s <%s>\nSubject: An Error occured generating Barcodes.\n\n    
    This is the error:
    
    %s
    
    %s
    """%('Thomas van Gurp', 'thomas@deenabio.com',opts, error)
    try:
       smtpObj = smtplib.SMTP('localhost')
       smtpObj.sendmail(sender, receivers[0], message)         
       print """A download link for the Barcodes was sent to you by email
       
       Thanks for your patience!"""
    except SMTPException:
       print """Error: unable to send email
       
       Please provide a valid email-adress"""

def write_log(opts, time):
    """Writes the failed / succesful attempt to the log file"""
    log_file = "../httpdocs/GBS-Barcodeusage_log.csv"
    time = time.localtime()
    with open(log_file, "a") as log:
        csvwriter_log = csv.writer(log,delimiter=';',quoting=csv.QUOTE_MINIMAL)
        csvwriter_log.writerow([time.tm_mday, time.tm_mon, time.tm_year, "%s:%s"%(time.tm_hour, time.tm_min),  \
        opts.name, opts.organization, opts.email, opts.enzyme,opts.number, opts.output ])
    
def run_script(opts):
    """runs barcodev2.py from the cgi script using subprocess"""
    args = ' '
    args += '--enzyme %s '%opts.enzyme
    args += '-n %s '%opts.number
    args += '--min %s '%opts.min
    args += '--max %s '%opts.max
    args += '--mono_nt %s '%opts.mono_nt
    args += '-a %s '%opts.ad_seq
    args += '--common_adapter_sequence %s '%opts.cm_ad_seq
    args += '--output %s '%opts.output
#    args += '--email %s'%opts.email
#    args += '--name %s '%opts.name
#    args += '--organization %s'%opts.organization
    cmd = 'python ./barcodev2.py %s'%args
#    print cmd
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell= True)
    stdout, stderr = proc.communicate()
#    print stdout, stderr
    if stderr:
        return stderr
#    print stdout, stderr
    proc.wait()
    return 0

def parse_options(form):
    """Parse Options Provided in command-line"""
    parser = OptionParser()
    parser.add_option("-e", "--enzyme",  metavar = "enzyme",  
                      action = "store", type="string", default= form['Enzyme'].value,
                      dest = "enzyme", help="Choose restriction enzyme")
    parser.add_option("-n", "--number",  metavar = "number", 
                      default= form['number'].value,
                      type = "int", action = "store", 
                      help = "Number of output Barcodes")
    parser.add_option("--min",  metavar = "min",  action = "store", 
                      default = form['min'].value, dest = "min", type = "int",
                      help = "Minimum Barcode length")
    parser.add_option("--max",  metavar = "max", default = form['max'].value,
                       type = "int", action = "store", dest = "max", 
                       help = "Maximum Barcode length")
    parser.add_option("-m", "--mono_nt",  metavar = "mono_nt", 
                       default = form['mono_nt'].value,
                       type = "int", action = "store", dest = "mono_nt", 
                       help = "Maximum number of mononucleotides in Barcode")
    parser.add_option("-a","--barcode_adapter_sequence",  metavar = "ad_seq",
                      action = "store", dest = "ad_seq",
                      default = form['ad_seq'].value,                      
                      help = "Choose alternative 5'-->3'sequence for barcode adapter")
    parser.add_option("--common_adapter_sequence",  metavar = "cm_ad_seq",
                      action = "store", dest = "cm_ad_seq",
                      default = form['cm_ad_seq'].value,                      
                      help = "Choose alternative 5'-->3'sequence for common adapter")
    parser.add_option("-o", "--output",  metavar = "output",  
                      action = "store", type="string",default = '../httpdocs/barcode_output/',  
                      dest = "output", help = "Output file")
    parser.add_option( "--email",  metavar = "email",  
                      action = "store", type="string",default = form['email'].value,  
                      dest = "email", help = "email address")
    parser.add_option( "--name",  metavar = "name",  
                      action = "store", type="string",default = form['name'].value,  
                      dest = "name", help = "username")
    parser.add_option( "--organization",  metavar = "organization",  
                      action = "store", type="string",default = form['organization'].value,  
                      dest = "organization", help = "organization")
    opts, args = parser.parse_args()
    return opts, args

opts,args = parse_options(form)
opts.output += str("%s-barcodes-%s-size-%s-%s-v2.0.csv"%(opts.number, opts.enzyme, opts.min, opts.max))
return_code = run_script(opts)
if return_code == 0:
    print """content-type: text/html\n\n

    <!doctype html public '-//W3C//DTD HTML 4.01//EN'
      'http://www.w3.org/TR/html4/strict.dtd'>
    <html>
      <head>
        <title>GBS-Barcode 2.0</title>
      </head>
      <H1>Processing complete!</H1>
      <p>please check your email for the download link</p>
      <p>Please note that the method with which barcodes are generated has been changed significantly.
      <br>
      The script is >100x faster and some bugs are fixed. 
      <br><br>
      If you have any questions regarding
      the new version please contact us at <a href='mailto:thomas@deenabio.com'>deenabio</a>"""
    email_success(opts)
    write_log(opts, time)
elif 'MemoryError' in return_code:
    print """content-type: text/html\n\n

    <!doctype html public '-//W3C//DTD HTML 4.01//EN'
      'http://www.w3.org/TR/html4/strict.dtd'>
    <html>
      <head>
        <title>GBS-Barcode 2.0</title>
      </head>
      <H1>Processing failed!</H1>
      <p>Due to a memory error the script unfortunately could not generate the
      specified number of barcodes. This issue is being adressed currently. In the meantime, please try a lower maximum length for your barcodes.</p>

      <br><br>
      If you have any questions regarding
      the new version please contact us at <a href='mailto:thomas@deenabio.com'>deenabio</a>"""
    write_log(opts, time)
    email_fail(opts, return_code)
elif 'cut in' in return_code:
    index = return_code.index('[') +1
    position =re.sub("\D", "", return_code[index:index+2])
    if 'common adapter' in return_code:
        adapter = opts.cm_ad_seq
    else:
        adapter = opts.ad_seq
    rb = Restriction.CommOnly
    ad_seq = Seq(adapter, IUPACAmbiguousDNA())
    analysis = Restriction.Analysis(rb, ad_seq)
    print """content-type: text/html\n\n
    <!doctype html public '-//W3C//DTD HTML 4.01//EN'
      'http://www.w3.org/TR/html4/strict.dtd'>
    <html>
      <head>
        <title>GBS-Barcode 2.0</title>
      </head>
      <H1>Processing failed!</H1>
      <p>
        Enzyme: <a href='http://rebase.neb.com/rebase/enz/%s.html'>%s</a> will cut (one of) the adapter sequences, see below.
       </p>
        It is <b>strongly</b> adviced to choose another enzyme.
      <br><br>
      If you have any questions regarding this
      please contact us at <a href='mailto:thomas@deenabio.com'>deenabio</a><PLAINTEXT>"""%(opts.enzyme, opts.enzyme)
    analysis.print_as('map')
    if 'common adapter' in return_code:
        print analysis.print_that(title = 'The following enzymes will cut the common adapter:\n\n')
    else:
        print analysis.print_that(title = 'The following enzymes will cut the barcoded adapter:\n\n')
    write_log(opts, time)
else:
    error_message =  """content-type: text/html\n\n

    <!doctype html public '-//W3C//DTD HTML 4.01//EN'
      'http://www.w3.org/TR/html4/strict.dtd'>
    <html>
      <head>
        <title>GBS-Barcode 2.0</title>
      </head>
      <H1>Processing failed!</H1>
      <p>An unknown error occured. The error will be investigated.
      <br><br>
      If you'd like to come back to this please contact us at <a href='mailto:thomas@deenabio.com'>deenabio</a>"""
    email_fail(opts, return_code)



