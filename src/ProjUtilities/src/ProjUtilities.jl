#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\src\ScoreDrivenERGM.jl\ProjUtilities.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\src\ScoreDrivenERGM.jl
#Created Date: Saturday April 24th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Saturday April 24th 2021 6:12:38 pm
#Modified By:  Domenico Di Gangi
#-----
#Description:
#-----
########

module ProjUtilities

using TableView
using Blink

viewtab(df) = body!(Window(), showtable(df))
export viewtab

end