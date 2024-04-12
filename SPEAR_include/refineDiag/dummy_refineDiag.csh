#!/bin/csh

echo ""
echo "  -- begin dummy_refineDiag.csh --  "
echo ""
#The mere existance of a refineDiag section in the xml pointing to any non-empty refineDiag
#script (as simple as doing a "echo" above)
#causes the history files to be unpacked by frepp in /ptmp/$USER/$ARCHIVE/$year.nc
#when the current year data lands on gfdl archive.
#
echo "  ---------- begin yearly analysis ----------  "
echo ""

echo "  ---------- end yearly analysis ----------  "

echo "  -- end   dummy_refineDiag.csh --  "
