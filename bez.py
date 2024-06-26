from scipy.special import binom
import matplotlib.pyplot as plt
import math
from geom import geom
from vector import vector
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
        
        # the first is the weight of the point
        
        # all beyond that define the spatial dimesions
        
        # a few simple curves would be
        
        # [[1,0],[1,1]] defining a line from 0 to 1 in 1 d space
        # [[1,0,0],[1,1,1]] defining a line from (0,0) to (1,1) in 2d space
        
        # and so on
        
        self.control_points=control_points
    
    def discretize2(self,order=1,stepsize=0.1,absolute=None):
        """full analytic representation of the curve keep orders lower 
        if you want performance"""
        
        #ah. shiet. can't decide on which points influence which
        #position without introducing another relative measure.
        
        #yes I can. I can use the original stepsize along the linear
        #elements
        #or more precisely, I can divide the length of the vector between
        #two control points with my stepsize and I get a counter of how 
        #many samples I should take from where.
        
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve
        c=0
        m=1
        
        divs=1/len(self.control_points)
        point_l=len(self.control_points)
        #influencing points= if c+order*divs<point_div_val<c+order*divs
        #eh. no.
        # rel_dist I want to evaluate
        # equal space distances of points derived from number of points
        # get point that's closest to rel_dist,
        # then get it's index and count in both directions for order steps
        # have those influence the point with the old order formula.
        
        
        dims=len(self.control_points[0])-1
        points=[]#self.control_points[0][2:]
        #this generates the running variable between 0 and 1
        
        #order=len(self.control_points)-1
        #ah ok, I'm tying order to length of my list.
        #higher list, higher order
        
        m=len(self.control_points)
        
        #I need the total length first,
        mysum=0
        v_lengths=[]
        c=0
        m=len(self.control_points)
        while c < m-1:
            p1x=self.control_points[c]
            p1=vector.Vector(*p1x,0)
            p2x=self.control_points[c+1]
            p2=vector.Vector(*p2x,0)
            pf=p2-p1
            
            
            
            c+=1
            li=pf.magnitude()
            mysum+=li
            v_lengths.append(li)
        
        c=0
        
        next_i=0
        
        while c < m:
           
            point=[0]*dims
            oc = 0
            
            i=0
            ns=0
            for x in v_lengths:
                ns+=x
                
                if ns>c*mysum:
                    break
                i+=1
                
            i=int(i)
            
            rel_point_c=i
            
            end=min(i+order+1,m-1)
            
            
            relevant_points=self.control_points[rel_point_c:end]
            if i==end:
                break
                
            #this is the sum from the formula
            while oc <= order and len(relevant_points)>order:
                #rint(oc,order)
                b=binom(order,oc)


                #shouldn't the position of the point
                #in my controlpointlist have an influence?
                #i.e. bigger order goes further outside?
                #of course, if I'm looking at the curve as a whole
                #I have to take the entire curve into account.
                
                #everything else needs a different approach.
                #
                #influence per point
                v=(1-c)**(order-oc)*(c)**oc
                #dimensionen
                
                dimc=0
                #these are the dimensions of the points
                
                #eh=int(c+oc)
                
                while dimc < dims:
                    
                    #self.control_points[eh]
                    dimv=b*v*relevant_points[oc][dimc+1]
                   
                    
                    point[dimc]+=dimv
                    dimc+=1
                    
                oc+=1
                
            
            if absolute:
                
                
                #ok erm how about reducing c by a fraction of the stepsize
                
                
                if len(points)>0:
                    d=distance(points[-1],point)
                    low = absolute - 0.1 * absolute
                    high = absolute + 0.1 * absolute
                    
                    if round(stepsize,5)==0:
                        
                        break
                    if  low < d < high:
                        points.append(point)
                        c+=stepsize
                        continue
                    
                    elif d < low :
                        c-=stepsize
                        
                        stepsize=stepsize*1.1
                        
                        if round(d)==0:
                            break
                            
                    elif d > high:
                        c-=stepsize
                                                
                        stepsize=stepsize*0.9
                        
                    c+=stepsize 
                else:
                    points.append(point)
                    
                    c+=stepsize
            else:
                points.append(point)
            #points.append(self.point(c))
                c+=stepsize
        return points
    def make_factor_index_d(self):
        """
        builds a factor from 0 to 1, that maps to the location
        of the point in the control point list.
        """
        factor_index_d = {}
        c_i = 0
        m_fs = len(self.control_points)
        while c_i < m_fs:
            factor = c_i/(m_fs-1)
            factor_index_d[factor] = c_i
            c_i += 1
        return factor_index_d
    
    def make_point(self, c, order):
        dims = len(self.control_points[0])
        point = [0] * dims
        oc = 0
        #this is the sum from the formula
        while oc <= order:
            b = binom(order,oc)
            
            v = (1-c)**(order-oc)*(c)**oc
            dimc = 0
            while dimc < dims:
                dimv = b * v * self.control_points[oc][dimc]
                point[dimc] += dimv
                dimc += 1
            oc += 1
        return point
    
    
    def get_other_point(self,factor_index_d,c):
        cp_factors = list(factor_index_d.keys())
        cp_factors.sort()
        c_list = 0
        m_list = len(cp_factors)-1
        while c_list < m_list:
            low = cp_factors[c_list]
            high = cp_factors[c_list+1]
            if low <= c < high:
                break
            c_list += 1
        
        l_p1 = list(self.control_points[c_list])
        l_p2 = list(self.control_points[c_list+1])
        if len(l_p1)<3:
            l_p1 += [0]
            l_p2 += [0]
        l_p1 = vector.Vector(*l_p1)
        l_p2 = vector.Vector(*l_p2)
        
        d_l = high-low
        d_c = c-low
        other_point = l_p2*(d_c/d_l)+l_p1*(1-(d_c/d_l))
        return other_point
    
    def bezier_to_linear(self,stepsize_factor=1,linear_factor=0, stepsize=None):
        """
        the idea behind this, is that I want a smooth-ish transition
        in a series of bezier curves that keeps the step size factor
        but moves closer to what the linear interpolation between
        the points looks like.
        """
        """full analytic representation of the curve keep orders lower 
        if you want performance"""
        
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve
        
        
        dims = len(self.control_points[0])
        if stepsize ==None:
            stepsize = 1/(len(self.control_points)*4)
        points = [] 
        
        order = len(self.control_points)-1
        # ah ok, I'm tying order to length of my list.
        # higher list, higher order
    
        factor_index_d = self.make_factor_index_d()
        
        # the curve goes from 0 to 1
        c = 0
        m = 1
        while c <= m:
            # make the smooth point
            point = self.make_point(c, order)
            # make the not smooth point
            other_point = self.get_other_point(factor_index_d,c)
            
            if len(point)<3:
                point=(*point,0)
            point = vector.Vector(*point)
            
            # interpolate between them
            inf1 = linear_factor * other_point
            inf2 = (1-linear_factor) * point
            actual_point = inf1 + inf2
            
            points.append(actual_point)
            
            c += stepsize
            c = round(c,3)
                
        return points
    
    def discretize(self,stepsize_factor=1,absolute=None):
        """full analytic representation of the curve keep orders lower 
        if you want performance"""
        
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve
        c = 0
        m = 1
        
        dims=len(self.control_points[0])-1
        stepsize=1/(len(self.control_points)*4)
        points=[]#self.control_points[0][2:]
        #this generates the running variable between 0 and 1
        
        order=len(self.control_points)-1
        #ah ok, I'm tying order to length of my list.
        #higher list, higher order
        
        while c <= m:
           
            point=[0]*dims
            oc = 0
            
            #this is the sum from the formula
            while oc <= order:
                b=binom(order,oc)
                
                v=(1-c)**(order-oc)*(c)**oc
                #dimensionen
                dimc=0
                #these are the dimensions of the points
                
                while dimc < dims:
                    dimv=b*v*self.control_points[oc][dimc+1]
                   
                    point[dimc]+=dimv
                    dimc+=1
                    
                oc+=1
                
            
            
            points.append(point)
            c+=stepsize
            c=round(c,3)
                
        return points

def test_plot_discretize2(c,stepsize,label="relative",color=None,order=1):
    x=[]
    y=[]
    for p in c.discretize2(order,stepsize):
        x.append(p[0])
        y.append(p[1])
    if color==None:
        plt.plot(x,y,label=label)
    else:
        plt.plot(x,y,label=label,color=color)
def test_plot(c,stepsize,label="relative",color=None):
    x=[]
    y=[]
    for p in c.discretize(stepsize):
        x.append(p[0])
        y.append(p[1])
    if color==None:
        plt.plot(x,y,label=label)
    else:
        plt.plot(x,y,label=label,color=color)

def make_normals(c,stepsize,label="relative"):
    x=[]
    y=[]
    points=c.discretize(stepsize)
    
    c=0
    m=len(points)
    normals={}
    while c < m-1:
        p1=points[c]
        p2=points[c+1]
        bdx=p2[0]-p1[0]
        bdy=p2[1]-p1[1]
        
        midx=p1[0]+bdx/2
        midy=p1[1]+bdy/2
        
        dx=bdx/2
        dy=bdy/2
        n=(-dy,dx)
        np1=(midx,midy)
        np2=(midx+n[0],midy+n[1])
        normals[c]=(np1,np2)
        c+=1
    return normals

def test_plot_normals(normals):
    keyl=list(normals.keys())
    keyl.sort()
    for key in keyl:
        tup=normals[key]
        np1,np2=tup
    
    
        plt.plot([np1[0],np2[0]],[np1[1],np2[1]],color="blue")
    #plt.plot(x,y,color="blue")    
    #plt.plot(x,y,label=label)
    
def test_plot_absolute(c,stepsize,absolute_distance):
    x=[]
    y=[]
    for p in c.discretize(stepsize,absolute_distance):
        x.append(p[0])
        y.append(p[1])
    plt.plot(x,y,label="absolute")



    
def test_simple(stepsize=0.1):
    ps2=[[1,1,0],[1,1,1],[1,2,1],[1,2,0],[1,3,0]]
    c2=Curve(ps2)
    ps,other=c2.discretize()
    x=[]
    y=[]
    for x in ps:
        x.append(p[0])
        y.append(p[1])
    plt.plot(x,y,label=label,color=color)

def test_absolute():
    ps=[[0,1,0,0],[0.5,1,0,1],[1,1,15,1]]
    c=Curve(ps)
    ste=0.1
    test_plot(c,ste)
    test_plot_absolute(c,ste,1)
    plt.xlim([0,20])
    plt.ylim([0,2])
    for x in ps:
        plt.plot([x[2]],[x[3]],"r+")
    plt.axis('equal')
    plt.legend()
    plt.show()

def make_water_lines():
    wls=[0.5,1,1.25,1.5,1.75]
    wls=[0.5,1.25,1.5]
    for wl in wls:
        ps=[[1,15,0],[1,12.5,wl],[1,10,wl],[1,3,wl],[1,0,wl*(3/5)]]
        yield wl,ps

def plot_waterlines():
    for tup in make_water_lines():
        ps=tup[1]
        wl=tup[0]
        #for correct orientation with x -> 
        ps.reverse()
        c=Curve(ps)
        ste=0.1
        test_plot(c,ste,str(wl))
        #test_plot_absolute(c,ste,1)
        normals=make_normals(c,ste)
        test_plot_normals(normals)
        plt.xlim([0,20])
        plt.ylim([0,2])
        for x in ps:
            plt.plot([x[1]],[x[2]],"r+")
            
    

def plot_spanten():
    ste=0.1
    #this roundness roughly 1m-2m from bow
    zs=[[1,0,0],[1,1,0.2],[1,1.2,2]]
    zs_curve=Curve(zs)
    rel_ps=zs_curve.discretize(stepsize=0.05)
    test_plot(zs_curve,ste)
    
    zs=[[1,0,0],[1,1,0.2],[1,1.2,2]]
    zs_curve=Curve(zs)
    rel_ps=zs_curve.discretize2(stepsize=0.05,order=1)
    test_plot(zs_curve,ste,color="red")
    
    #plt.plot([0,1],[0,0.2])
    #plt.plot([1,1.2],[0.2,2])
    
    
def plot_contour():
    
    #startpunkt wäre die intersection der spline mit n(0,1,0) 
    #bei dem y wert den ich angebe.
    #oder einer der werte die ich beim plotten benutze die y>wert haben.
    #dann vllt. noch interpolieren, dann wäre ich relativ nahe dran.
    
    xs=[0,1,2,3,4,5]
    
    ste=0.1
    zs=[1,0,0],[1,2,0],[1,4,1]
    zs=[[1,0,1],[1,3,0],[1,5,0],[1,10,0],[1,12,0],[1,14,1]]
    plt.plot([14,15],[1,2.5],color="blue")
    plt.plot([0,0.5],[1,2.5],color="blue")
    plt.plot([0.5,15],[2.5,2.5],color="blue")
    zs_curve=Curve(zs)
    rel_ps=zs_curve.discretize(stepsize=0.05)
    test_plot(zs_curve,ste)
    
    
    return
    
    
def lower_order_long_test():
    ste=0.1
    #zs=[[1,0,1],[1,3,0],[1,5,0],[1,10,0],[1,12,0],[1,14,1]]
    
    zs=[[1,15,0],[1,12.5,1],[1,10,1],[1,3,1],[1,0,1*(3/5)]]
    zs_curve=Curve(zs)
    #original order=len(zs)-1
    test_plot_discretize2(zs_curve,ste,order=4)
    plt.show()
    
def boat_stuff():
    
    
    #waterlines
    plot_waterlines()
    #plot_contour()
    #plot_spanten()
    #ok, spanten?
    plt.axis('equal')
    plt.legend()
    plt.show()
    
    #anyway, what do I do with this?
    
    # I spawn voronoi cells along the splines
    # then check angles and areas
    # hmmm.


def easier_test():
    s=[[1,15,0],[1,12.5,1],[1,10,1],[1,3,1],[1,0,1*(3/5)]]
    
    s.reverse()
    a=s[0:4]
    b=s[3:]
    spline_sequence=[a,b]
    for s in spline_sequence:
        C=Curve(s)
        test_plot(C,0.1,color="blue",label="easy?")
    s=[[1,15,0],[1,12.5,1],[1,10,1],[1,3,1],[1,0,1*(3/5)]]
    s.reverse()
    C=Curve(s)
    test_plot(C,0.1,color="red",label="easy?")
    plt.axis('equal')
    plt.legend()
    plt.show()
    #plt.show()

def interp_f(x,mul=5):
    inp=x*mul
    return 1-abs(inp)

def interp():
    #ok so.
    #bezier curves always incompass the entire controlpoint list
    #I don't want that. And it would be cool if the curve passed through
    #the control points 
    #I sometimes just want a bezier - like curve that's based on the
    #control point environment, but I also don't want to use splines.
    #I think.
    
    #these will still not go through the definig points
    
    ps2=[[1,1,0],[1,1,1],[1,2,1],[1,2,0],[1,3,0]]
    
    r=1/len(ps2)
    
    #now to define my interpolation function.
    #it has a range somewhere between 0 and 1, probably a lot smaller
    
    #ok I have one function defined, now I need to select the applying control
    #points for any given 0<c<1
    
    c=0
    m=len(ps2)
    while c <= 1:
        
        ri=interp_f(c,m)
        
        c+=r
        c=round(c,3)
    
    a=1

def test_bez_lin():
    cps=[(0,0),(0,1),(1,1)]#,(1,0),(2,0)]
    C=Curve(cps)
    fs=[0,0.25,0.5,0.75,1]
    for f in fs:
        
        disc=C.bezier_to_linear(linear_factor=f)
        x=[x[0] for x in disc]
        y=[x[1] for x in disc]
        plt.plot(x,y,color=(f,0,0))
    plt.axis("equal")
    plt.show()

def test_geom():
    points = [(0,0),(0,-1),(1,-0.9),(1,-2.1),(-1,-1.9),(-1,-3.1)]
    curve = Curve(points)
    these = curve.bezier_to_linear()
    geoml = []
    c = 0
    while c < len(these)-1:
        #print(these[c])
        p1 = vector.Vector(these[c][0],these[c][1],0)
        p2 = vector.Vector(these[c+1][0],these[c+1][1],0)
        l = geom.Line.from_two_points(p1,p2)
        geoml.append(l)
        c += 1
    fl = []
    for x in geoml:
        fl.append(x.as_svg())
    
    view_box_d = geom.make_view_box_d(geoml)
    geom.main_svg(fl,"beztest.svg",view_box_d=view_box_d)


class EasySpline:
    def __init__(self,controlpoints):
        """AAAAAAAA THIS ISN'T EXACT YET."""
        # ok so the original ideas would be that we're calculating
        # coefficients
        # Since we can calculate normal and normal vectors in 2d and 3d,
        # we don't actually need to do that.
        self.controlpoints = controlpoints
    
    def interpolate(self,steps=5,order=1):
        """idk if this works for more than 1"""
        # so.what we are doing, specifically.
        # well, the order specifices continuity types
        # 0, the lines are continues in the 0th derivative. Which just means they meet.
        # 1 the lines are continues in the 1st derivative, which is slope.
        # this containts the clue how to solve this, because we can take
        
        
        # calculate the difference vectors between them.
        
        self.lists_of_order_vectors = [self.controlpoints]
        order_c = 0
        while order_c < order:
            lower_list = self.lists_of_order_vectors[order_c]
            
            diffs = build_diffs(lower_list)
            new_list = build_derivatives(diffs)
            self.lists_of_order_vectors.append(new_list)
            order_c+=1
            
        output_points = []
        point_c = 0
        m = len(self.controlpoints)
        while point_c < m-1:
            print("")
            p1 = self.controlpoints[point_c]
            p2 = self.controlpoints[point_c+1]
            print(self.lists_of_order_vectors[order])
            diff1 = self.lists_of_order_vectors[order][point_c]
            diff2 = self.lists_of_order_vectors[order][point_c+1]
            inter_c = 0
            while inter_c < steps:
                
                # how about the point definitions is
                
                # p1 + diff and that's decreasing the closer I get
                # and the other ons i p2-diff and that's increasing the closer I get.
                # I think that's what I'm already doing?
                # but it still can't be done on the next first step.
                
                inf1 = (1-(inter_c/(steps-1)))
                inf2 = (inter_c/(steps-1))
                
                pr = p1 * inf1 + p2 * inf2
                div1 = diff1 * inf1 * inf2
                print(diff1,div1)
                pr += div1 
                div2 = diff2 * inf2 * inf1
                pr -= div2
                
                
                print(inter_c,p1,p2,pr)
                
                output_points.append(pr)
                inter_c += 1
            point_c += 1 
        return output_points

def build_derivatives(diffs):
    new_list = []
    point_c = 0
    point_m = len(diffs)
    while point_c < point_m:
        if point_c == 0:
            new_list.append(diffs[0])
        #elif point_c == point_m-1:
            #new_list.append(diffs[-1])
        else:
            diff1 = diffs[point_c-1]
            diff2 = diffs[point_c]
            deriv = (diff1 + diff2)/2
            new_list.append(deriv)
        point_c += 1
    
    new_list.append(diffs[-1])
    print("new list")
    print(new_list)
    
    return new_list

def build_diffs(lower_list):
    diffs = []
    point_c = 0
    point_m = len(lower_list)
    while point_c < point_m-1:
        p1 = lower_list[point_c]
        p2 = lower_list[point_c+1]
        diff = (p2-p1)
        diffs.append(diff)
        point_c += 1
    return diffs

def test_easy_spline():
    
    p1 = vector.Vector(0,0,0)
    p2 = vector.Vector(1,1,0)
    p3 = vector.Vector(2,0,0)
    p4 = vector.Vector(3,1,0)
    p5 = vector.Vector(4,0,0)
    p6 = vector.Vector(5,-1,0)
    
    cps = [p1,p2,p3,p4,p5]#,p6]
    ES = EasySpline(cps)
    
    geoml = []
    for x in ES.controlpoints:
        geoml.append(geom.circle(x*3,radius=0.1))
    output = ES.interpolate(order=1)
    
    c=0
    m= len(ES.lists_of_order_vectors[1])
    while c < m:
        p=ES.controlpoints[c]*3
        v=ES.lists_of_order_vectors[1][c]*3
        line=geom.Line.from_two_points(p-0.2*v,p+0.2*v)
        line.style ="stroke:rgb(255,0,0);"
        geoml.append(line)
        c+=1
    #for x in ES.lists_of_order_vectors[1]
    
    these = output
    c = 0
    while c < len(these)-1:
        #print(these[c])
        p1 = vector.Vector(these[c][0],these[c][1],0)*3
        p2 = vector.Vector(these[c+1][0],these[c+1][1],0)*3
        l = geom.Line.from_two_points(p1,p2)
        geoml.append(l)
        c += 1
    fl = []
    for x in geoml:
        fl.append(x.as_svg())
    
    view_box_d = geom.make_view_box_d(geoml)
    geom.main_svg(fl,"EZ_Spline_test.svg",view_box_d=view_box_d)

if __name__=="__main__":
    #easier_test()
    #
    #r=test_simple()
    #plt.axis('equal')
    #plt.show()
    #lower_order_long_test()
    #boat_stuff()
    
    #interp()
    
    #test_bez_lin()
    #test_geom()
    test_easy_spline()
