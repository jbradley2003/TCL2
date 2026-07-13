import numpy as np

#
#
### 2X2 SPIN OPERATORS - spin 1/2
#
#

sx = np.array(([[0,1/2],[1/2,0]]),dtype=np.complex128)
sy = np.array(([[0,-1/2],[1/2,0]]),dtype=np.complex128)*1j
sz = np.array(([[1/2,0],[0,-1/2]]),dtype=np.complex128)
Id = np.identity(2)

splus = np.array([[0,1],[0,0]],dtype=np.complex128)
sminus = np.array([[0,0],[1,0]],dtype=np.complex128)

#
#
### 3x3 SPIN OPERATORS - spin 1
#
#

sx3 = np.array(([[0,1,0],[1,0,1],[0,1,0]]),dtype=np.complex128)*(1/np.sqrt(2))
sy3 = np.array(([[0,-1,0],[1,0,-1],[0,1,0]]),dtype=np.complex128)*(1j/np.sqrt(2))
sz3 = np.array(([[1,0,0],[0,0,0],[0,0,-1]]),dtype=np.complex128)
Id3 = np.identity(3)

s3plus = np.array([[0,1,0],[0,0,1],[0,0,0]],dtype=np.complex128)*np.sqrt(2)
s3minus = np.array([[0,0,0],[1,0,0],[0,1,0]],dtype=np.complex128)*np.sqrt(2)


#
#
### 4x4 SPIN OPERATORS - spin 3/2, etc.
#
#

sx4 = np.array(([[0,np.sqrt(3),0,0],[np.sqrt(3),0,2,0],[0,2,0,np.sqrt(3)],[0,0,np.sqrt(3),0]]),dtype=np.complex128)*(1/2)
sy4 = np.array(([[0,-np.sqrt(3),0,0],[np.sqrt(3),0,-2,0],[0,2,0,-np.sqrt(3)],[0,0,np.sqrt(3),0]]),dtype=np.complex128)*(1j/2)
sz4 = np.array(([[3/2,0,0,0],[0,1/2,0,0],[0,0,-1/2,0],[0,0,0,-3/2]]),dtype=np.complex128)
Id4 = np.identity(4)

s4plus = np.array([[0,np.sqrt(3),0,0],[0,0,2,0],[0,0,0,np.sqrt(3)],[0,0,0,0]],dtype=np.complex128)
s4minus = np.array([[0,0,0,0],[np.sqrt(3),0,0,0],[0,2,0,0],[0,0,np.sqrt(3),0]],dtype=np.complex128)

#
#
### 5x5 SPIN OPERATORS - spin 2
#
#

sx5 = np.array(([[0,2,0,0,0],[2,0,np.sqrt(6),0,0],[0,np.sqrt(6),0,np.sqrt(6),0],[0,0,np.sqrt(6),0,2],[0,0,0,2,0]]),dtype=np.complex128)*(1/2)
sy5 = np.array(([[0,-2,0,0,0],[2,0,-np.sqrt(6),0,0],[0,np.sqrt(6),0,-np.sqrt(6),0],[0,0,np.sqrt(6),0,-2],[0,0,0,2,0]]),dtype=np.complex128)*(1j/2)
sz5 = np.array(([[2,0,0,0,0],[0,1,0,0,0],[0,0,0,0,0],[0,0,0,-1,0],[0,0,0,0,-2]]),dtype=np.complex128)
Id5 = np.identity(5)

s5plus = sx5 + 1j*sy5
s5minus = sx5 - 1j*sy5

#
#
### 6x6 SPIN OPERATORS - spin 5/2
#
#

sx6 = np.array(([[0,np.sqrt(5),0,0,0,0],[np.sqrt(5),0,np.sqrt(8),0,0,0],[0,np.sqrt(8),0,np.sqrt(9),0,0],[0,0,np.sqrt(9),0,np.sqrt(8),0],[0,0,0,np.sqrt(8),0,np.sqrt(5)],[0,0,0,0,np.sqrt(5),0]]),dtype=np.complex128)*(1/2)
sy6 = np.array(([[0,-np.sqrt(5),0,0,0,0],[np.sqrt(5),0,-np.sqrt(8),0,0,0],[0,np.sqrt(8),0,-np.sqrt(9),0,0],[0,0,np.sqrt(9),0,-np.sqrt(8),0],[0,0,0,np.sqrt(8),0,-np.sqrt(5)],[0,0,0,0,np.sqrt(5),0]]),dtype=np.complex128)*(1j/2)
sz6 = np.array(([[5/2,0,0,0,0,0],[0,3/2,0,0,0,0],[0,0,1/2,0,0,0],[0,0,0,-1/2,0,0],[0,0,0,0,-3/2,0],[0,0,0,0,0,-5/2]]),dtype=np.complex128)
Id6 = np.identity(6)

s6plus = sx6 + 1j*sy6
s6minus = sx6 - 1j*sy6

#
#
### 7x7 SPIN OPERATORS - spin 3
#
#

sx7 = np.array([[0,np.sqrt(3/2),0,0,0,0,0],
               [np.sqrt(3/2),0,np.sqrt(5/2),0,0,0,0],
               [0,np.sqrt(5/2),0,np.sqrt(3),0,0,0],
               [0,0,np.sqrt(3),0,np.sqrt(3),0,0],
               [0,0,0,np.sqrt(3),0,np.sqrt(5/2),0],
               [0,0,0,0,np.sqrt(5/2),0,np.sqrt(3/2)],
               [0,0,0,0,0,np.sqrt(3/2),0]],dtype=np.complex128)

sy7 = np.array([[0,-np.sqrt(3/2),0,0,0,0,0],
               [np.sqrt(3/2),0,-np.sqrt(5/2),0,0,0,0],
               [0,np.sqrt(5/2),0,-np.sqrt(3),0,0,0],
               [0,0,np.sqrt(3),0,-np.sqrt(3),0,0],
               [0,0,0,np.sqrt(3),0,-np.sqrt(5/2),0],
               [0,0,0,0,np.sqrt(5/2),0,-np.sqrt(3/2)],
               [0,0,0,0,0,np.sqrt(3/2),0]],dtype=np.complex128)*1j

sz7 = np.array([[3,0,0,0,0,0,0],
               [0,2,0,0,0,0,0],
               [0,0,1,0,0,0,0],
              [0,0,0,0,0,0,0],
               [0,0,0,0,-1,0,0],
               [0,0,0,0,0,-2,0],
               [0,0,0,0,0,0,-3]],dtype=np.complex128)

Id7 = np.identity(7)

s7plus = sx7 + 1j*sy7
s7minus = sx7 - 1j*sy7

#
#
### 8x8 SPIN OPERATORS - spin 7/2
#
#

sx8 = np.array(([[0,np.sqrt(7),0,0,0,0,0,0],[np.sqrt(7),0,np.sqrt(12),0,0,0,0,0],[0,np.sqrt(12),0,np.sqrt(15),0,0,0,0],[0,0,np.sqrt(15),0,np.sqrt(16),0,0,0],[0,0,0,np.sqrt(16),0,np.sqrt(15),0,0],[0,0,0,0,np.sqrt(15),0,np.sqrt(12),0],[0,0,0,0,0,np.sqrt(12),0,np.sqrt(7)],[0,0,0,0,0,0,np.sqrt(7),0]]),dtype=np.complex128)*(1/2)
sy8 = np.array(([[0,-np.sqrt(7),0,0,0,0,0,0],[np.sqrt(7),0,-np.sqrt(12),0,0,0,0,0],[0,np.sqrt(12),0,-np.sqrt(15),0,0,0,0],[0,0,np.sqrt(15),0,-np.sqrt(16),0,0,0],[0,0,0,np.sqrt(16),0,-np.sqrt(15),0,0],[0,0,0,0,np.sqrt(15),0,-np.sqrt(12),0],[0,0,0,0,0,np.sqrt(12),0,-np.sqrt(7)],[0,0,0,0,0,0,np.sqrt(7),0]]),dtype=np.complex128)*(1j/2)
sz8 = np.array(([[7/2,0,0,0,0,0,0,0],[0,5/2,0,0,0,0,0,0],[0,0,3/2,0,0,0,0,0],[0,0,0,1/2,0,0,0,0],[0,0,0,0,-1/2,0,0,0],[0,0,0,0,0,-3/2,0,0],[0,0,0,0,0,0,-5/2,0],[0,0,0,0,0,0,0,-7/2]]),dtype=np.complex128)
Id8 = np.identity(8)

s8plus = sx8 + 1j*sy8
s8minus = sx8 - 1j*sy8

spin_operators = {0.5:{'sx':sx,'sy':sy,'sz':sz,'Id':Id,'s+':splus,'s-':sminus},
                  1.0:{'sx':sx3,'sy':sy3,'sz':sz3,'Id':Id3,'s+':s3plus,'s-':s3minus},
                 1.5:{'sx':sx4,'sy':sy4,'sz':sz4,'Id':Id4,'s+':s4plus,'s-':s4minus},
                  2.0:{'sx':sx5,'sy':sy5,'sz':sz5,'Id':Id5,'s+':s5plus,'s-':s5minus},
                  2.5:{'sx':sx6,'sy':sy6,'sz':sz6,'Id':Id6,'s+':s6plus,'s-':s6minus},
                  3.0:{'sx':sx7,'sy':sy7,'sz':sz7,'Id':Id7,'s+':s7plus,'s-':s7minus},
                  3.5:{'sx':sx8,'sy':sy8,'sz':sz8,'Id':Id8,'s+':s8plus,'s-':s8minus}
                 }


#
#
### CONSTANTS
#
#

gyro_dict = {'e':-17608.59705,'H':26.75221824,'C':6.72828532,'V':7.0492,'N':1.93297,'Cu':7.1118,'P':10.829,'Cl':2.62401,'D':4.106645934475427, 'Cl':3.077,'Mn':6.641,'B':2.874,'Br':6.719,'F':25.176} # modify this as needed - in units of rad ms^-1 G^-1
spin_type_dict = {'e':0.5,'H':0.5,'C':0.5,'V':3.5,'N':1.0,'Cu':1.5,'P':0.5,'Cl':2,'Br':1.5,'D':1.0,'Mn':2.5,'B':3.0,'F':0.5,'1/2':0.5,'1':1.0,'3/2':1.5,'2':2.0,'5/2':2.5,'3':3.0,'7/2':3.5}
