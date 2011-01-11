# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   rspec -I ../bioruby/lib/ spec/gffdb_spec.rb 
#
# Copyright (C) 2010 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"


require 'bio-gff3'

include Bio::GFFbrowser

describe FastLineParser, "GFF3 Fast line parser" do
  include Bio::GFFbrowser::FastLineParser

  it "should parse attributes" do 
    parse_attributes_fast("id=1").should == { "id"=>"1" }
    parse_attributes_fast("id=1;parent=45").should == { "id"=>"1", "parent" => "45" }
  end
  it "should parse records" do 
    parse_line_fast("ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1").should ==  ["ctg123", ".", "CDS", 1201, 1500, 0.0, "+", 0, {"ID"=>"cds00001", "Parent"=>"mRNA00001", "Name"=>"edenprotein.1"}]
  end
end


