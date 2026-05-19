from mission_subsystems import Propulsion, AOCS, CaptureSystem, EPS, CDH, Structures, TCS, TTC

class SpaceCraft:
    def __init__(self, sc_type, m_wet_t=None, contingency=0.3):
        self.sc_type = sc_type
        self.contingency = contingency

        if self.sc_type == "ms":
            self.m_wet_t = m_wet_t

        # dry mass independent subsystems
        self.aocs = AOCS(self.sc_type)
        self.m_aocs = self.aocs.mass()

        self.capture_system = CaptureSystem(self.sc_type)
        self.m_capture_system = self.capture_system.mass()

        self.eps = EPS(self.sc_type)
        self.m_eps = self.eps.mass()

        self.cdh = CDH()
        self.m_cdh = self.cdh.mass()

        self.dry_mass()
        self.prop_mass()
        
    
    def dry_mass(self):
        '''Estimate the dry mass of the spacecraft'''
        m_dry = self.m_aocs + self.m_capture_system + self.m_eps + self.m_cdh  # dry mass independent systems
        m_dry_prev = 0
        while abs(m_dry - m_dry_prev) > 1:
            m_dry_prev = m_dry

            if self.sc_type == "ms":
                propulsion = Propulsion(m_dry, self.sc_type, m_wet_t=self.m_wet_t) 
            else:
                propulsion = Propulsion(m_dry, self.sc_type) 
            structures = Structures(m_dry)
            tcs = TCS(m_dry)
            ttc = TTC(m_dry)

            m_propulsion = propulsion.mass()
            m_structures = structures.mass()
            m_tcs = tcs.mass()
            m_ttc = ttc.mass()

            m_dry = (self.m_aocs + self.m_cdh + self.m_capture_system + self.m_eps) + m_propulsion + m_structures + m_tcs + m_ttc

        if self.sc_type == "ms":
            self.propulsion = Propulsion(m_dry, self.sc_type, m_wet_t=self.m_wet_t) 
        else:  # for the tug
            self.propulsion = Propulsion(m_dry, self.sc_type)
        self.structures = Structures(m_dry)
        self.tcs = TCS(m_dry)
        self.ttc = TTC(m_dry)

        self.m_propulsion = self.propulsion.mass()
        self.m_structures = self.structures.mass()
        self.m_tcs = self.tcs.mass()
        self.m_ttc = self.ttc.mass()
        self.m_dry = m_dry

        return m_dry 
    
    def prop_mass(self):
        m_prop_item = self.propulsion.propellant_mass()
        m_prop = m_prop_item[0].mass_kg

        self.m_prop = m_prop
        return m_prop

    def mass_breakdown(self):
        return {
            "mass_aocs": self.m_aocs,
            "mass_capture_system": self.m_capture_system,
            "mass_eps": self.m_eps,
            "mass_cdh": self.m_cdh,
            "mass_propulsion": self.m_prop,
            "mass_structures": self.m_structures,
            "mass_tcs": self.m_tcs,
            "mass_ttc": self.m_ttc
        }
