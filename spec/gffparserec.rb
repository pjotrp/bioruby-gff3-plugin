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

  it "should parse attributes" do 
  end
end


