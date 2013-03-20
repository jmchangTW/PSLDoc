#!/usr/bin/env perl

open (OUT_F, ">config.h");
print OUT_F "#ifndef _CONFIG_H\n";
print OUT_F "#define _CONFIG_H\n\n";

print "Process PSLDoc.cfg...\n";
open (F, "PSLDoc.cfg");
while(!eof(F))
{
	$line = <F>;

	if((!($line =~ m/^#.+/))&&(length($line) > 2))
	{
		my @arr = split(/\s/, $line);
		if($arr[1] eq "=")
		{		
			output_define($arr[0],$arr[2]);
			if(($arr[0] eq PSIBLAST_DATABASE_PATH) and (length($arr[2]) > 1))
			{
				 `export BLASTDB=$arr[2]`;
				 print "specify the path of blast search databse = $arr[2]\n";
			}
		}
	}
}
close (F);
print "Sucess\n";
print OUT_F "#endif /* _CONFIG_H */\n";
print "Output configure information to \"config.h\" for \"make\"\n";

sub output_define()
{	
	my($key, $value) = @_;
	print OUT_F "#ifndef $key\n";
	print OUT_F "#define $key \"$value\"\n";
	print OUT_F "#endif\n\n";
}

