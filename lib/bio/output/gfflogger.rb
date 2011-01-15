
module Bio
  module GFFbrowser

    module Helpers

      module Logger
        include Bio::Log

        def debug str, id=''
          log = LoggerPlus['bio-gff3']
          log.debug str+" <#{id}>"
        end

        def info str, id=''
          log = LoggerPlus['bio-gff3']
          log.info str+" <#{id}>"
        end

        def warn str, id=''
          log = LoggerPlus['bio-gff3']
          log.warn str+" <#{id}>"
        end

        def error str, id=''
          log = LoggerPlus['bio-gff3']
          log.error_(str+" <#{id}>",:act => FailOnError.new)

        end

        def log_sys_info msg
          log = LoggerPlus['bio-gff3']
          rmem = `ps -o rss= -p #{Process.pid}`.to_i
          vmem = `ps -o vsz= -p #{Process.pid}`.to_i
          if rmem or vmem
            log.info7 "Memory used #{msg} RAM #{rmem/1024}M, VMEM #{vmem/1024}M"
          end
        end
      end
    end
  end
end
