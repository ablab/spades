#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import sendgrid


def send_message(user_list, subject, text):
    sg = sendgrid.SendGridClient("spades", "5phGnoqf2n")
    message = sendgrid.Mail()

    message.add_to("spades.support@bioinf.spbau.ru")
    for address in user_list:
        message.add_bcc(user_list)

    message.set_from("spades.support@bioinf.spbau.ru")
    message.set_subject(subject)
    message.set_html(text)

    sg.send(message)


def save_user_list(file_name, email_list):
    out_file = open(file_name, 'w')
    email_set = set(email_list)
    
    for email in sorted(email_set):
        out_file.write(email + '\n')
    out_file.close()


if len(sys.argv) != 3:
    print('Mail sender usage: ' + sys.argv[0] + ' <file with emails> <file with messge>')
    print('Emails should be stored one per line')
    print('File with message should contain email subject in the first line and message in HTML format afterwards.')
    sys.exit(0)

user_list_file = open(sys.argv[1])
user_list = []
for address in user_list_file:
    if address.find('@') != -1 and address.find('.') != -1:
        user_list.append(address.strip())

user_set = set(user_list)
if len(user_set) < len(user_list):
    save = raw_input('Warning! Found ' + str(len(user_list) - len(user_set)) + ' duplicate enties. Do you wish to save filtered list? (y/n) ')
    if (save.upper().startswith('Y')):
        out_file_name = 'filtered_list.txt'
        custom_name = raw_input('Enter file name (default name is ' + out_file_name + '): ')
        if custom_name != '':
            out_file_name = custom_name
        save_user_list(out_file_name, user_set)

message_file = open(sys.argv[2])
subject = message_file.readline().strip()
message = ''
for line in message_file:
    message += line

send = raw_input('Are you sure you want to send the message with subject "' + subject + '" to ' + str(len(user_set)) + ' users? (yes/no) ')
if (send.strip().upper() == 'YES'):
    send_message(user_set, subject, message)
    print("Message sent")
else:
    print("Aborted")


