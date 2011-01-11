
module Bio
  module GFFbrowser

    module Helpers

      class Section < Range
        attr_reader :rec
        def initialize rec
          super(rec.start,rec.end)
          @rec = rec
        end
        def intersection(other)
          raise ArgumentError, 'value must be a Range' unless other.kind_of?(Range)
          min, max = first, exclude_end? ? max : last
          other_min, other_max = other.first, other.exclude_end? ? other.max : other.last
          new_min = self === other_min ? other_min : other === min ? min : nil
          new_max = self === other_max ? other_max : other === max ? max : nil
          new_min && new_max ? new_min..new_max : nil
        end
        alias_method :&, :intersection
        def <=> other
          first <=> other.first
        end
      end

      module Sections
        # Return list of sorted Sections
        def Sections::sort rec
          sections = []
          rec.each do | section |
            sections.push Section.new(section)
          end
          sections.sort
        end
      end

    end
  end
end

