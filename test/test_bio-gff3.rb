$: << '.'
do_create = if ARGV[0] == '-c' or ARGV[0] == '--create'
              ARGV.shift
            end
  
require 'test/unit'
require 'regressiontest'

RegressionTest.create(do_create)

class Gff3Test < Test::Unit::TestCase

  rootpath = File.join(File.dirname(__FILE__),'..')
  BIN = rootpath + '/bin/gff3-fetch'
  DAT = rootpath + '/test/data'
  def test_cache
    assert_equal(true,single_run("mRNA #{DAT}/gff/test-ext-fasta.fa #{DAT}/gff/test-ext-fasta.gff3",'test_ext_gff3'))
    assert_equal(true,single_run("CDS #{DAT}/gff/test.gff3",'test_gff3'))
  end

  def test_nocache
    assert_equal(true,single_run("mRNA --cache none #{DAT}/gff/test-ext-fasta.fa #{DAT}/gff/test-ext-fasta.gff3",this_method+'_ext_gff3'))
    assert_equal(true,single_run("CDS #{DAT}/gff/test.gff3",this_method+'_gff3'))
  end

  private
   def this_method
     caller[0] =~ /`([^']*)'/ and $1
   end

end

def single_run opts, name
  cmd = "#{BIN} --logger stdout #{opts}"
  # p cmd
  RegressionTest.test `#{cmd}`,name,"#{DAT}/regression"
end
