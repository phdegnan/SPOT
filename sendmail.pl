#!/usr/bin/perl
##
## 	Patrick Degnan
##	sendmail.pl
##	Send notice of job completion 
##

unless(@ARGV){die "\nUsage: $0 job_id jane.doe\@email.com originating_program.pl\n\n";}

$job=$ARGV[0];
$email=$ARGV[1];
$program=$ARGV[2];

$send_to = "To: $email\n"; 
$sendmail = "/usr/sbin/sendmail -t"; 	
#$reply_to = "Reply-to: your\@email.edu\n"; 
$subject = "Subject: $program job $job complete\n";

$content="
Dear User,

$program job $job has completed. Login to retrieve your results.\n

Thanks,
Admin

";	
	
	

open(SENDMAIL, "|$sendmail") or die "Cannot open $sendmail: $!";

print SENDMAIL $reply_to; 
print SENDMAIL $subject; 
print SENDMAIL $send_to; 
print SENDMAIL $cc;
print SENDMAIL "Content-type: text/plain\n\n"; 
print SENDMAIL $content; 

print "$reply_to $subject $send_to $cc Content-type: text/plain\n\n $content";


close(SENDMAIL); 

