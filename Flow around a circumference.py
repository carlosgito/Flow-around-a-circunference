#Import of the libraries
import math as math
import matplotlib.pyplot as plt
import numpy as np

#Definition of the constants
n = 200     #number of rows and columns (our grid always will be a square)
V = 37.5    #velocity m/s
L = 200     #length (m)
dn = L/n    #length of each row and each column
rho = 1.225
radius = 50

#Definition of the vectors and matrixes
vectorXaxis = np.arange(0,L+1, dn)
vectorYaxis = np.arange(0,L+1, dn)
array = [[0 for x in range(n+1)] for y in range(n+1)]   #stream function
array2 = [[0 for x in range(n+1)] for y in range(n+1)]  #ball
arrayV = [[0 for x in range(n+1)] for y in range(n+1)]  #y-direction velocity
arrayU = [[0 for x in range(n+1)] for y in range(n+1)]  #x-direction velocity
arrayMv = [[0 for x in range(n+1)] for y in range(n+1)] #velocity module
arrayP = [[0 for x in range(n+1)] for y in range(n+1)]  #pressure

#Numeric definition of the circle. It will check if each point
#of the matrix is inside or a boundary of the circle:
N=0
while N<len(array):
    M=0
    while M<len(array):
        r = math.sqrt((((M-(len(array)/2))**2)+(N-(len(array)/2))**2))
        #All the points outside the circle get a 0 value in the array2 matrix
        #All the points inside the circle get a 1 value in the array2 matrix
        if r<radius:
            array2[N][M] = 1
        #All the points at the boundary of the circle get a 2 value in the array2 matrix
        if r==radius:
            array2[N][M] = 2
        M += 1
    N += 1

#Boundary conditions. We will use the ones described in the project's first part
k=0
while k <= n:
    array[0][k] = 0.0
    k += 1
k = 0
while k <= n:
    array[n][k] = V*L
    k += 1
k = 0
while k <= n:
    array[k][n] = V*vectorYaxis[k]
    k += 1
k=0
while k <= n:
    array[k][0] = V*vectorYaxis[k]
    k += 1

#Calculation of the stream function values.
#We will use 5000 iterations to make it more accurate
arraynew = array*1
k = 1
j = 1
iterations = 1

while iterations <= 5000:
    k = 1
    while k < n:
        j = 1
        while j < n:
            formula = ((arraynew[j + 1][k] + arraynew[j - 1][k]) + (arraynew[j][k + 1] + (arraynew[j][k - 1])))/4
            #The program makes a difference between the points inside the circle,
            #in the boundary or outside
             #Outside
            if array2[j][k] == 0:
                array[j][k] = round(formula, 4)
             #Inside (it theorically is =0, but we will use the same as in the boundary
             #to increase the efficiency in each iteration.
            if array2[j][k]==1:
                formula = V*L/2
                array[j][k] = round(formula, 4)
             #Boundary
            if array2[j][k] == 2:
                formula = V*L/2
                array[j][k] = round(formula, 4)
            j += 1
        k += 1
    arraynew = array*1
    iterations = iterations + 1

#X-component of the velocity:
k=1
while k < n:
    j = 1
    while j<n:
        deltaPsiX=array[j][k]-array[j-1][k]
        incrementoDX = j * dn - (j - 1) * dn
        u=deltaPsiX/incrementoDX
        #The normal velocity must be 0 in the boundary of the ball (as boundary condition)
        #and inside the ball
        if array2[j][k] == 2 or array2[j][k] == 1:
            arrayU[j][k]=0
        else:
            arrayU[j][k]=round(u,3)
        j+=1
    k+=1

#Y-component of the velocity:
j=1
while j<n:
    k = 1
    while k < n:
        deltaPsiY= array[j][k-1] - array[j][k]
        incrementoDY = k * dn - (k - 1) * dn
        v=deltaPsiY/incrementoDY
        #Velocity inside and in the boundary of the ball must be 0
        if array2[j][k] == 2 or array2[j][k] == 1:
            arrayV[j][k]=0
        else:
            arrayV[j][k]=round(v,3)
        k+=1
    j+=1

#Velocity module
j=1
while j<n:
    k = 1
    while k < n:
        vm = math.sqrt((arrayU[j][k]) ** 2 + (arrayV[j][k]) ** 2)
        arrayMv[j][k]=round(vm,4)
        k+=1
    j+=1


#Boundary Conditions of the velocity and the pressure
k=0
while k <= n:
    arrayU[0][k] = V
    arrayV[0][k] = 0
    arrayP[0][k] = 101350
    k += 1
k = 0
while k <= n:
    arrayU[k][n] = V
    arrayV[k][n] = 0
    arrayP[k][n] = 101350
    k += 1
k = 0
while k <= n:
    arrayU[n][k] = V
    arrayV[n][k] = 0
    arrayP[n][k] = 101350
    k += 1
k=0
while k <= n:
    arrayU[k][0] = V
    arrayV[k][0] = 0
    arrayP[k][0] = 101350
    k += 1

#Pressure
j=1
vo = ((math.sqrt((arrayU[0][1])**2+(arrayV[0][1])**2))**2)
while j<n:
    k = 1
    while k < n:
        P = arrayP[0][1] + (rho/2)*(vo - (math.sqrt((arrayU[j][k])**2+(arrayV[j][k])**2))**2)
        arrayP[j][k]=round(P,4)
        k+=1
    j+=1

#Plotting of the different graphics
plt.title("Stream function in m^2/s")
plt.contourf(array, 20)
plt.grid(color='white', linewidth=0.5)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
circle1=plt.Circle((len(array)/2,len(array)/2),radius, facecolor='black')
plt.gcf().gca().add_artist(circle1)
plt.gca().set_aspect('equal', adjustable='box')
plt.axis(xlim=(0,200), ylim=(0,200))
plt.colorbar()
plt.show()

plt.title("Velocity in the x axis in m/s")
plt.contourf(arrayU, 20)
plt.grid(color='white', linewidth=0.5)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
circle1=plt.Circle((len(array)/2,len(array)/2),radius, facecolor='black')
plt.gcf().gca().add_artist(circle1)
plt.gca().set_aspect('equal', adjustable='box')
plt.colorbar()
plt.show()

plt.title("Velocity in the y axis in m/s")
plt.contourf(arrayV, 20)
plt.grid(color='white', linewidth=0.5)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
circle1=plt.Circle((len(array)/2,len(array)/2),radius, facecolor='black')
plt.gcf().gca().add_artist(circle1)
plt.gca().set_aspect('equal', adjustable='box')
plt.colorbar()
plt.show()

plt.title("Velocity module m/s")
plt.contourf(arrayMv, 20)
plt.grid(color='white', linewidth=0.5)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
circle1=plt.Circle((len(array)/2,len(array)/2),radius, facecolor='black')
plt.gcf().gca().add_artist(circle1)
plt.gca().set_aspect('equal', adjustable='box')
plt.colorbar()
plt.show()

plt.title("Pressure Field in Pa")
plt.contourf(arrayP, 20)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.grid(color='white', linewidth=0.5)
circle1=plt.Circle((len(array)/2,len(array)/2),radius, facecolor='black')
plt.gcf().gca().add_artist(circle1)
plt.gca().set_aspect('equal', adjustable='box')
plt.colorbar()
plt.show()
