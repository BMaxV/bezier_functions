from scipy.special import binom
import math

from vector import vector


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


class EasySpline:
    def __init__(self,controlpoints):
        """AAAAAAAA THIS ISN'T EXACT YET."""
        # ok so the original ideas would be that we're calculating
        # coefficients
        # Since we can calculate normal and normal vectors in 2d and 3d,
        # we don't actually need to do that.
        self.controlpoints = controlpoints
    
    def interpolate(self,steps=10,order=1):
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
            p1 = self.controlpoints[point_c]
            p2 = self.controlpoints[point_c+1]
            diff1 = self.lists_of_order_vectors[order][point_c]
            diff2 = self.lists_of_order_vectors[order][point_c+1]
            inter_c = 0
            while inter_c < steps:
                
                regular = (inter_c/(steps-1))
                regular_inv = 1-(inter_c/(steps-1))
                
                p1_inf = p1 + diff1 * regular * regular_inv
                p2_inf = p2 - diff2 * regular * regular_inv
                
                sine_fac_down = (math.cos(regular*math.pi)+1)/2
                sine_fac_up = (math.sin(regular*math.pi - math.pi/2)+1)/2
                
                pr = p1_inf * sine_fac_down + p2_inf * sine_fac_up 
                
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
