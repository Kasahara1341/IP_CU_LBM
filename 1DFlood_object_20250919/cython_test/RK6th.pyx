class Runge_Kutta_6th(object):
    def __init__(self) -> None:
        self.stage = 0
        self.hold = None
        self.increments = []
        
    def update_stage(self):
        if self.stage < 5:
            self.stage += 1
        else:
            self.stage = 0
    
    def set_hold(self,value):
        self.hold = value
    
    def update_stage0_variables(self,Element,dt):
        k1 = Element.calc_increment()
        uppdated_depth = self.hold+0.2*k1*dt
        Element.set_depth(uppdated_depth)
        self.increments.append(k1)
        self.update_stage()
    
    def update_stage1_variables(self,Element,dt):
        k2 = Element.calc_increment()
        uppdated_depth = self.hold+(3/40)*self.increments[0]*dt+(9/40)*k2*dt
        Element.set_depth(uppdated_depth)
        self.increments.append(k2)
        self.update_stage()

    # ここから
    def update_stage2_variables(self,Element,dt):
        k3 = Element.calc_increment()
        uppdated_depth = self.hold+(0.3)*self.increments[0]*dt-0.9*self.increments[1]*dt+1.2*k3*dt
        Element.set_depth(uppdated_depth)
        self.increments.append(k3)
        self.update_stage()

    def update_stage3_variables(self,Element,dt):
        k4 = Element.calc_increment()
        uppdated_depth = self.hold-(11/54)*self.increments[0]*dt+2.5*self.increments[1]*dt-(70/27)*self.increments[2]*dt+(35/27)*k4*dt
        Element.set_depth(uppdated_depth)
        self.increments.append(k4)
        self.update_stage()

    def update_stage4_variables(self,Element,dt):
        k5 = Element.calc_increment()
        uppdated_depth = self.hold+(1631/55296)*self.increments[0]*dt+(175/512)*self.increments[1]*dt+(575/13824)*self.increments[2]*dt+(44275/110592)*self.increments[3]*dt+(253/4096)*k5*dt
        Element.set_depth(uppdated_depth)
        self.increments.append(k5)
        self.update_stage()

    def update_stage5_variables(self,Element,dt):
        k6 = Element.calc_increment()
        self.increments.append(k6)
        self.update_stage()
    

    def update_depth(self,Element,dt):
        if self.stage == 0:
            self.hold = Element.get_variable_depth()
            self.update_stage0_variables(Element,dt)

        elif self.stage==1:
            self.update_stage1_variables(Element,dt)
        
        elif self.stage==2:
            self.update_stage2_variables(Element,dt)

        elif self.stage==3:
            self.update_stage3_variables(Element,dt)

        elif self.stage==4:
            self.update_stage4_variables(Element,dt)

        else:
            self.update_stage5_variables(Element,dt)
            uppdated_depth = self.hold + ((37/378)*self.increments[0]+0*self.increments[1]+(250/621)*self.increments[2]+(125/594)*self.increments[3]+(0)*self.increments[4]+(512/1771)*self.increments[5])*dt
            Element.set_depth(uppdated_depth)
            self.increments = []