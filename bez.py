from scipy.special import binom
import matplotlib.pyplot as plt
import math
#ok... patches?
#simple curve works

def distance(p1,p2):
    if len(p1)==len(p2) and len(p1) < 3:
        p1=[p1[0],p1[1],0]
        p2=[p2[0],p2[1],0]
    """the euclidian 3d distance between two points"""
    d=math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)
    return d


class Curve:
    def __init__(self,control_points):
        #control_points is a list of coordinates.
        
        # the first is the position on the curve,
        # the second is the weight of the point
        
        # all beyond that define the spatial dimesions
        
        # a few simple curves would be
        
        # [[0,1,0],[1,1,1]] defining a line from 0 to 1 in 1 d space
        # [[0,1,0,0],[1,1,1,1]] defining a line from (0,0) to (1,1) in 2d space
        
        # and so on
        
        self.control_points=control_points
    
    
    def discretize(self,stepsize=0.1,absolute=None):
        """full analytic representation of the curve keep orders lower 
        if you want performance"""
        
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve
        c=0
        m=1
        points=[]
        #this generates the running variable between 0 and 1
        while c < m:
            print("c",c)
            #print("looping through my steps",c)
            dims=len(self.control_points[0])-2
            point=[0]*dims
            oc = 0
            order=len(self.control_points)-1
            #this is the sum from the formula
            while oc <= order:
                #print("order",oc)
                b=binom(order,oc)
                #print("binom",b)
                v=(1-c)**(order-oc)*c**oc
                #print("einfluss",v)
                #dimensionen
                
                dimc=0
                #these are the dimensions of the points
                
                #print("cp",self.control_points[oc])
                while dimc < dims:
                    dimv=b*v*self.control_points[oc][dimc+2]
                   
                    point[dimc]+=dimv
                    dimc+=1
                
                oc+=1
            #print("\n point",point,"\n")
            if absolute:
                print("absolute",absolute)
                if len(points)>0:
                    d=distance(points[-1],point)
                    low = absolute - 0.1 * absolute
                    high = absolute + 0.1 * absolute
                    print("d,low,high",d,low,high)
                    if round(stepsize,5)==0:
                        print("breaking")
                        break
                    if  low < d < high:
                        points.append(point)
                        print("good stepsize")
                        c+=stepsize
                    
                    elif d < low :
                        print(stepsize)
                        c=c*1.1
                        stepsize=stepsize*1.1
                        print(stepsize)
                        if round(d)==0:
                            break
                        print(d)
                        print(low)
                        print("bad stepsize, too low, increasing")
                        
                    elif d > high:
                        print(stepsize)
                        c=c*0.9
                        stepsize=stepsize*0.9
                        print(stepsize)
                        print("bad stepsize, too high, lowering")
                        
                    else:
                        print("woot")
                else:
                    points.append(point)
                    
                    c+=stepsize
            else:
                points.append(point)
            #points.append(self.point(c))
                c+=stepsize
        return points

def test_plot(c,stepsize):
    x=[]
    y=[]
    for p in c.discretize(stepsize):
        #print(round(p[0],3),round(p[1],3))
        x.append(p[0])
        y.append(p[1])
    plt.plot(x,y,label="relative")
def test_plot_absolute(c,stepsize,absolute_distance):
    x=[]
    y=[]
    for p in c.discretize(stepsize,absolute_distance):
        #print(round(p[0],3),round(p[1],3))
        x.append(p[0])
        y.append(p[1])
    plt.plot(x,y,label="absolute")



def test():
    #ps1=[[0,1,0],[1,1,1]]
    #c1=Curve(ps1)
    #test_plot(c1)
    test_simple()
    
def test_simple(stepsize=0.1):
    ps2=[[0,1,0,0],[0.5,1,0,1],[1,1,1,1]]
    c2=Curve(ps2)
    test_plot(c2,stepsize)
    plt.show()

def test_absolute():
    ps=[[0,1,0,0],[0.5,1,0,1],[1,1,15,1]]
    c=Curve(ps)
    ste=0.15
    test_plot(c,ste)
    test_plot_absolute(c,ste,0.1)
    plt.legend()
    plt.show()


if __name__=="__main__":
    #test()
    test_absolute()
    
