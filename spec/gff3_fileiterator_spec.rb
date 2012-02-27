# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   rspec -I ../bioruby/lib/ spec/gff3_fileiterator_spec.rb 
#
# Copyright (C) 2010 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"

require 'bio-gff3'

TEST1='test/data/gff/test.gff3'
TEST2='test/data/gff/standard.gff3'

describe Bio::GFF::GFF3::FileIterator, "iterates a GFF3 file" do


  it "should parse a file and yield records" do 
    iter = Bio::GFF::GFF3::FileIterator.new(TEST1)
    iter.each_rec do | fpos, line |
      rec = Bio::GFF::GFF3::FastParserFileRecord.new(fpos, line)
      rec.io_seek.should == 51
      break
    end
  end

  it "should handle embedded FASTA records" do
    iter = Bio::GFF::GFF3::FileIterator.new(TEST1)
    last = nil
    iter.each_rec do | fpos, line |
      rec = Bio::GFF::GFF3::FastParserFileRecord.new(fpos, line)
      last = rec
    end
    last.io_seek.should == 3342
    firstid = []
    iter.each_sequence do | id, seq |
      firstid << id
    end
    firstid.sort.should == ['test01','test02']
  end

end


