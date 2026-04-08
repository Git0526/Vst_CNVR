#!/usr/bin/perl
use strict;
use Data::Dumper;
die "Usage: perl $0 <CNV table> <case ID> <control ID> <output> <plus>\n" if(@ARGV <4);
=cut 
Note:
1. Vst =(Vtotal–[Vpop1×Npop1+Vpop2×N pop2]/Ntotal)/Vtotal

Vtotal：pop1和pop2的方差

V：Variance of pop
N：Number of pop 

Vst was calculated using the equation: Vst =(Vtotal–[Vpop1×Npop1+Vpop2×Npop2]/Ntotal)/Vtotal, 
where Vtotal is the total variance, 
Vpop is the variance for each respective population, 
Npop is the sample size for each respective population, 
and Ntotal is the total sample size. 
Subsequently, the CNVRs with the top five Vst values were taken as candidate CNVRs, 
and functional enrichment analysis of these regions was performed. 

2. Vs = [Vpop1×Npop1+Vpop2×N pop2]/Ntotal


=cut 
open (IN1, "$ARGV[0]") or die "CNV number table required!\n";
open (IN2, "$ARGV[1]") or die "case ID required!\n";
open (IN3, "$ARGV[2]") or die "control ID required!\n";
open (OUT, ">$ARGV[3]") or die "permission denied!\n";
my $pos = $ARGV[4] if (@ARGV ==5);
my $total_sample = 0;
my %case;
while (<IN2>){
	chomp;
	$total_sample++;
	$case{$_} = 1;
}
close IN2;
#print Dumper(\%case);

my %control;
while (<IN3>){
	chomp;
	$total_sample++;
	$control{$_} = 1;
}
close IN3;
print STDERR "case control file loading finished!\n";
print STDERR "total sample number $total_sample\n";
my $header = <IN1>; chomp $header;
my @tmp = split/\s+/, $header;
#print OUT "$header\tVST\n";
my $index = -1; 
my (@case, @control,@case_t,@control_t);
foreach my $ID (@tmp){
	$index++;
	#next if $index > $total_sample+8;
	if ($ID  =~ /(\S+)_type/){
		if (exists $case{$1}){push @case_t, $index}
		if (exists $control{$1}){push @control_t, $index}
	}else{
		if (exists $case{$ID}){push @case, $index}
		if (exists $control{$ID}){push @control, $index}
	}
}
print STDERR "case individual number: ".scalar @case. "\ncontrol individual number: ".scalar @control."\n";
print OUT join("\t","#CHROM","POS","ID","REF","ALT","QUAL",@tmp[@control_t],@tmp[@case_t],@tmp[@control],@tmp[@case],"VST")."\n";
while (<IN1>){
	chomp;
	my @tmp_copynumber = split/\s+/;
	my @case_copynumber = @tmp_copynumber[@case];               ### depth_case
	my @control_copynumber = @tmp_copynumber[@control];         ### dep  control
	my @total_copynumber = @tmp_copynumber[@case, @control];    ### total
	my $total_variance = &var(\@total_copynumber);    
	my $case_variance = &var(\@case_copynumber);
	my $control_variance = &var(\@control_copynumber);
	my $weight_Vs = ($case_variance*($#case+1)+$control_variance*($#control+1))/($#case+$#control+2);
	my $VST = $total_variance > 0 ? sprintf("%.4f",($total_variance-$weight_Vs)/$total_variance) : 0;
	if (defined $pos){
		print OUT join("\t",@tmp_copynumber[0..5],@tmp_copynumber[@control_t],@tmp_copynumber[@case_t],@tmp_copynumber[@control],@tmp_copynumber[@case],$VST)."\n" if ($VST >=0);
	}else{
		print OUT join("\t",@tmp_copynumber[0..5],@tmp_copynumber[@control_t],@tmp_copynumber[@case_t],@tmp_copynumber[@control],@tmp_copynumber[@case],$VST)."\n";
	}
}
close OUT;

###方差 :   variance
sub var{
	my ($data) = @_;
	my $sum;
	for my $value (@{$data}) {
		$sum += $value;
	}
	my $mean = @{$data} > 0 ? $sum/@{$data} : 0;
	my $sqtotal = 0;
	foreach (@{$data}) {
		$sqtotal += ($_ - $mean) ** 2
	}
	my $var = @{$data} > 1 ? $sqtotal / (scalar @{$data} - 1) : 0;
	return $var;
}
