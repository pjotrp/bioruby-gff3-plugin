# To recreate the regression files:
#
#   ruby -Itest test/test_bio-gff3.rb --create

$: << '.'
do_create = if ARGV[0] == '-c' or ARGV[0] == '--create'
              ARGV.shift
            end
  
require 'test/unit'
require 'regressiontest'
require 'regressiontest2'

RegressionTest2.create(do_create)

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

  def test_latest_wormbase
    opts = "CDS #{DAT}/gff/m_hapla.WS232.genomic.part.fa #{DAT}/gff/m_hapla.WS232.annotations.part.gff3"
    arg1 = this_method + '_ext_gff3'

    bin = File.expand_path(BIN)
    cmd = "#{bin} --logger stdout #{opts}"
    assert_equal(true,RegressionTest::CliExec::exec(cmd,arg1))
  end

  private

  def this_method
    caller[0] =~ /`([^']*)'/ and $1
  end

end

def single_run opts, name
  bin = File.expand_path(BIN)
  cmd = "#{bin} --logger stdout #{opts}"
  if false
    print "Skipping ", cmd, "!\n"
  else
    text = `#{cmd}`.split(/\n/).delete_if { | s | s =~ /Memory/ }.join("\n")

    RegressionTest2.test text,name,"#{DAT}/regression"
  end
  true
end
