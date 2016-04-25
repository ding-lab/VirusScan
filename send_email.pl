#!/usr/bin/perl
use strict;
use warnings;
my ( $dir, $email ) = @ARGV;
my $sendmail = "/usr/sbin/sendmail -t"; 
my $reply_to = "Reply-to: $email\n"; 
my $subject = "Subject: data processing finished\n"; 
my $content = "The $dir data processing has finished.\n"; 
my $send_to = "To: $email\n"; 
open(SENDMAIL, "|$sendmail") or die "Cannot open $sendmail: $!"; 
print SENDMAIL $reply_to; 
print SENDMAIL $subject; 
print SENDMAIL $send_to; 
print SENDMAIL "Content-type: text/plain\n\n"; 
print SENDMAIL $content; 
close SENDMAIL;
exit;
