# To recreate the regression files:
#
#   ruby -Itest test/test_bio-gff3.rb --create

$: << '.'
do_create = if ARGV[0] == '-c' or ARGV[0] == '--create'
              ARGV.shift
            end
  
require 'test/unit'
require 'regressiontest'

RegressionTest.create(do_create)

rootpath = File.join(File.dirname(__FILE__),'..')
DAT = rootpath + '/test/data'
BIN = rootpath + '/bin/gff3-fetch'

class Gff3Test < Test::Unit::TestCase

  def test_cache
    assert_equal(true,single_run("mRNA #{DAT}/gff/test-ext-fasta.fa #{DAT}/gff/test-ext-fasta.gff3",'test_ext_gff3'))
    assert_equal(true,single_run("CDS #{DAT}/gff/test.gff3",'test_gff3'))
  end

  def test_nocache
    assert_equal(true,single_run("mRNA --cache none #{DAT}/gff/test-ext-fasta.fa #{DAT}/gff/test-ext-fasta.gff3",this_method+'_ext_gff3'))
    assert_equal(true,single_run("CDS --cache none #{DAT}/gff/test.gff3",this_method+'_gff3'))
  end

  def test_lrucache
    assert_equal(true,single_run("mRNA --cache lru #{DAT}/gff/test-ext-fasta.fa #{DAT}/gff/test-ext-fasta.gff3",this_method+'_ext_gff3'))
    assert_equal(true,single_run("CDS --cache lru #{DAT}/gff/test.gff3",this_method+'_gff3'))
  end

  private
   def this_method
     caller[0] =~ /`([^']*)'/ and $1
   end

end

def single_run opts, name
  bin = File.expand_path(BIN)
  cmd = "#{bin} --logger stdout #{opts}"
  p cmd
  # text = `#{cmd}`.split(/\n/).delete_if { | s | s =~ /Memory/ }.join("\n")

  # RegressionTest.test text,name,"#{DAT}/regression"
end
