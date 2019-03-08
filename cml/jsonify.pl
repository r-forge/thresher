#! perl -w

use strict;
use warnings;

my $source = shift or die "You must provide the name of a source file.";
my $target = "$source.json";
$target =~ s/.txt//;

print STDERR "$target\n";

open(SRC, "<$source") or die "Unable to open '$source': $!\n";
open(TGT, ">$target") or die "Unable to create '$target': $!\n";
print TGT "\{\n \"producer\": \"CytoGPS\",\n \"version\": \"1.0\",\n \"source\": \"$source\",\n"
    ." \"results\": [\n";

# use a manual finite state machine
my $state = "empty";
my $id = '';
my $msg = '';
my $clone = '';
my $indic = '';
my $LGF = '';
my $punct = '';

my $linenum = 0;
while (my $line = <SRC>) {
    ++$linenum;
    chomp($line);
    next unless $line; # skip blank lines
    if ($line =~ /^Event/ || $line =~ /^line/ || $line =~/^Maybe/) {
	$state = 'gunk';
	next;
    }
    if ($state eq "gunk") { # eat all lines until a  new entry starts
	next unless ($line =~ /^# (.*)$/); # New entry always starts with "# identifier"
	startGathering($1); # set state to 'starting' and define $id
	next;
    }
    if ($state eq 'empty') {
	if ($line =~ /^# (.*)$/) { # New entry always starts with "# identifier"
	    startGathering($1); # set state to 'starting' and define $id
	} else {
	    print STDERR "FAIL # $id $msg: $state\n";
	}
	next;
    }
    if ($state eq "starting") {
	if ($line =~ /^Karyotype (\d*)/) {
	    $clone = $1;
	    $state = "indicator";
	} else {
	    $msg .= "$line\n";
	}
	next;
    }
    if ($state eq "indicator") {
	if ($line =~ /\[(\d*)\]/) { # that's what we expect to see
	    $indic = $1;
	    $state = "LGF";
	} else {
	    die("Bad format (indicator, $linenum)\n");
	}
	next;
    }
    if ($state eq "LGF") {
	if ($line =~ /^\[(0|1), /) {
	    $LGF = $line;
	    $msg = "Passed" unless($msg);
	    # write the result
	    print TGT "$punct\{\n  \"id\": \"$id\",\n"
		."  \"message\": \"$msg\",\n"
		."  \"clone\": $clone,\n"
		."  \"indicator\": $indic,\n"
		."  \"LGF\": $LGF\n\}";
	    if ($indic > 0) {
		$LGF = 0;
		--$indic;
	    } else {
		# clear most variables
		$state = "tentative";
		$punct = ",\n";
		$LGF = $indic = $clone = $msg = '';
	    }
	} else {
	    die("Bad format (LGF; $linenum)\n");
	}
	next;
    }
    if ($state eq "tentative") {
	if ($line =~ /^Karyotype (\d*)/) {
	    $clone = $1;
	    $state = "indicator";
	} elsif ($line =~ /^# (.*)$/) { # New entry always starts with "# identifier"
	    startGathering($1); # set state to 'starting' and define $id
	} else {
	    die("Unexpected (tentative, $linenum\n");
	}
	next;
    }
    die "Unknown ($state, $linenum)): $line\n";
}

print TGT " ]\n\}\n";
close(TGT);
close(SRC);


sub startGathering {
    my $name = shift;
    $state = "starting";
    $id = $name;
}

exit;
__END__
