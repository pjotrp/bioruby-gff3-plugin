require 'rubygems'
require 'bundler'
begin
  Bundler.setup(:default, :development)
rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end
require 'rake'

require 'jeweler'
Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
  gem.name = "bio-gff3"
  gem.homepage = "https://github.com/pjotrp/bioruby-gff3-plugin"
  gem.license = "MIT"
  gem.summary = %Q{GFF3 parser for big data}
  gem.description = %Q{GFF3 (genome browser) information and digest mRNA and CDS sequences.
Options for low memory use and caching of records.
Support for external FASTA files.
}
  gem.email = "pjotr.prins@thebird.nl"
  gem.authors = ["Pjotr Prins"]
  # Include your dependencies below. Runtime dependencies are required when using your gem,
  # and development dependencies are only needed for development (ie running rake tasks, tests, etc)
  # gem.add_runtime_dependency 'bio', '>= 1.4.1'
  # gem.add_runtime_dependency 'log4r', '> 1.1.6'
  # gem.add_runtime_dependency 'bio-logger', '>= 0.8.0'
end
Jeweler::RubygemsDotOrgTasks.new

require 'rake/testtask'
Rake::TestTask.new(:test) do |test|
  test.libs << 'lib' << 'test'
  # test.pattern = 'test/**/test_*.rb'  # breaks in 1.9.3
  test.test_files = Dir.glob("test/**/test_*.rb")
  test.verbose = true
  Kernel.system('rspec spec/*.rb')
end

#require 'spec/rake/spectask'
#Spec::Rake::SpecTask.new(:spec) do |t|
#  t.spec_files = Dir.glob('spec/**/*_spec.rb')
#  t.spec_opts << '--format specdoc'
#  t.warning = true
#  t.rcov = true
#end


# require 'rcov/rcovtask'
# Rcov::RcovTask.new do |test|
#   test.libs << 'test'
#   test.pattern = 'test/**/test_*.rb'
#   test.verbose = true
# end

task :default => :test

# require 'rake/rdoctask'
# Rake::RDocTask.new do |rdoc|
#   version = File.exist?('VERSION') ? File.read('VERSION') : ""

#  rdoc.rdoc_dir = 'rdoc'
#   rdoc.title = "bio-gff3 #{version}"
#  rdoc.rdoc_files.include('README*')
#  rdoc.rdoc_files.include('lib/**/*.rb')
# end
