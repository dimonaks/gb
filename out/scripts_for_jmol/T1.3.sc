load "gb4/out/T1.3.xyz"
#set unitcell 1.4    #[line-width-or-type]
#Turns on or off the unit cell for crystal structures, and #determines its line style and line width (as a decimal #number, in Angstroms).

set boundbox 2
#Turns on or off a wire-frame box that contains the model, #and determines the line style and line width (as a decimal #number, in Angstroms) of that box.
background white
#rotate z 90
color [20,120,250]

select Ti1,Ti2,Ti3,Ti4,Ti5,Ti6,Ti7,Ti8,Ti9,Ti10,Ti11,Ti12,Ti13,Ti14,Ti15,Ti16,Ti17,Ti18,Ti19,Ti20,Ti21,Ti22,Ti23,Ti24,Ti25,Ti26,Ti27,Ti28,Ti29,Ti30,Ti31,Ti32,Ti33,Ti34,Ti35,Ti36,Ti69,Ti70,Ti71,Ti72,Ti73,Ti74,Ti75,Ti76,Ti77,Ti78,Ti79,Ti80,Ti81,Ti82,Ti83,Ti84,Ti85,Ti86,Ti87,Ti88,Ti89,Ti90,Ti91,Ti92,Ti93,Ti94,Ti95,Ti96,Ti97,Ti98,Ti99,Ti100,Ti101,Ti102,Ti103,Ti104
color red
select all
cpk 250
wireframe 0.3
#set axes 3


select * /2
translateSelected 0 -12 0
rotateSelected X 90
wireframe 0.1
cpk 150
model 0
zoom 60








#May usefull
#connect 6 #allow to connect with bond atoms on legthes less then 6 angstroms
#define mygroup within(5.0[FS4]102)
#select mygroup
#restrict
#show orientation

