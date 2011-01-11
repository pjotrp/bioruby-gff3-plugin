
module Bio
  module GFFbrowser

    module Helpers

      module Validate
        def validate_mrnas 
          return if not @options[:validate]
          # validate gene/container/component seqname is shared
          @mrnalist.validate_seqname
          @mrnalist.validate_shared_parent
        end

        def validate_cdss 
          return if not @options[:validate]
          @cdslist.validate_seqname
          # validate CDS sections do not overlap
          @cdslist.validate_nonoverlapping
          # validate sections share the parent
          @cdslist.validate_shared_parent
          # display unhandled features
        end
      end # Validate

      # Helper class for counting IDs
      class Counter < Hash
        def add id
          self[id] = 0 if self[id] == nil
          self[id] += 1
        end
      end

      # Helper class for storing linked records based on a shared ID
      class LinkedRecs < Hash
        include Error
        def add id, rec
          info "Adding #{rec.feature_type} <#{id}>"
          self[id] = [] if self[id] == nil
          self[id] << rec
        end

        # Validate all lists belong to the same container/component
        def validate_seqname
          each do | id, rec |
            seqname = rec.first.seqname
            rec.each do | section |
              raise "Non-matching seqname #{section.seqname} in #{seqname}" if section.seqname != seqname
            end
          end
        end
     
        # Validate all lists share the same parent (if available). First checks
        # for Parent attribute, next for mRNA attribute
        def validate_shared_parent
          each do | id, rec |
            parent = rec.first.get_attribute('Parent')
            if parent
              rec.each do | section |
                _parent = section.get_attribute('Parent')
                raise "Non-matching parent #{_parent} and #{parent} in #{id}" if _parent != parent
              end
            end
            parent = rec.first.get_attribute('mRNA')
            if parent
              rec.each do | section |
                _parent = section.get_attribute('mRNA')
                raise "Non-matching parent #{_parent} and #{parent} in #{id}" if _parent != parent
              end
            end
          end
        end
     
        # walk all (CDS) lists for every container/component and 
        # validate they do not overlap
        def validate_nonoverlapping
          each do | id, rec |
            sections = Sections::sort(rec)
            sections.each_with_index do | check, i |
              neighbour = sections[i+1]
              if neighbour and check.intersection(neighbour)
                warn "Overlapping sections for ",id
              end
            end
          end
        end
      end
    end
  end
end
