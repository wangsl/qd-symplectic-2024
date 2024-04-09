#!/bin/env perl

# $Id$

require 5.002;
use Getopt::Std;

$a = '$(O)/';

while(<>) {
    s/\/usr\S*\.h//g; 
    s# /usr/include/\S*##g;
    s# /usr/lib/\S*##g;
    #s/\/share\/apps\S*\.hpp//g; 
    #s/\/share\/apps\S*\.h//g; 
    s#/share/apps/\S*##g;
    #s/\w*\.cu//g;
    #s/\w*\.[Cc]//g;
    if (/\\$/) {
	s/\s*\\$//;
	chop;
    }
    s/ +/ /g;
    s/(\S+\.o)/$a$1/g;
    print;
} 

#$allfiles = join " ", grep !/rcsid.C/, glob( "*.[hcCfF]" );
#print "rcsid.C: $allfiles\n";
#print "\t perl rcsid.pl > rcsid.C\n";
