from subsystems import Propulsion, AOCS, Structures, TCS, CaptureSystem, EPS, TTC, CDH

class SpaceCraft:
    """
    SpaceCraft digital twin class
    """
    # initial phase: 35%, preliminary phase: 30%, detailed phase: 20%
    def __init__(self, n_targets, dv_list, m_debris_list, Isp, sc_type, n_redundancy, capture_type, mission_type, contingency=0.3):  # parameters needed for all the subsystem classes here
        self.n_targets = n_targets
        self.dv_list = dv_list
        self.m_debris_list = m_debris_list
        self.Isp = Isp
        self.sc_type = sc_type
        self.n_redundancy = n_redundancy
        self.capture_type = capture_type
        self.mission_type = mission_type
        self.contingency = contingency

        self.aocs = AOCS()
        self.cdh = CDH(self.n_redundancy)  
        self.capture_system = CaptureSystem(self.capture_type)
        self.eps = EPS(self.mission_type)

        self.m_aocs = self.aocs.mass()
        self.m_cdh = self.cdh.mass()
        self.m_capture_system = self.capture_system.mass()
        self.m_eps = self.eps.mass()
        

    def mass(self):
        """
        Compute the digitaltwin of the spacecraft

        """
        m_dry = self.m_aocs + self.m_cdh + self.m_capture_system + self.m_eps  # dry mass independent systems
        m_dry_prev = 0
        print(m_dry)
        while abs(m_dry - m_dry_prev) > 1:
            m_dry_prev = m_dry

            propulsion = Propulsion(m_dry, self.n_targets, self.dv_list, self.m_debris_list, self.Isp) 
            structures = Structures(m_dry)
            tcs = TCS(m_dry)
            ttc = TTC(m_dry)

            m_propulsion = propulsion.mass()
            m_structures = structures.mass()
            m_tcs = tcs.mass()
            m_ttc = ttc.mass()

            m_dry = (self.m_aocs + self.m_cdh + self.m_capture_system + self.m_eps) + m_propulsion + m_structures + m_tcs + m_ttc
            print(m_dry)

        self.propulsion = Propulsion(m_dry, self.n_targets, self.dv_list, self.m_debris_list, self.Isp) 
        self.structures = Structures(m_dry)
        self.tcs = TCS(m_dry)
        self.ttc = TTC(m_dry)

        self.m_propulsion = self.propulsion.mass()
        self.m_structures = self.structures.mass()
        self.m_tcs = self.tcs.mass()
        self.m_ttc = self.ttc.mass()
        self.m_dry = m_dry

        return m_dry * (1 + self.contingency)

    def mass_breakdown(self):
        self.m_aocs_percent = self.m_aocs/self.m_dry
        self.m_cdh_percent = self.m_cdh/self.m_dry
        self.m_capture_system_percent = self.m_capture_system/self.m_dry
        self.m_eps_percent = self.m_eps/self.m_dry
        self.m_propulsion_percent = self.m_propulsion/self.m_dry
        self.m_structures_percent = self.m_structures/self.m_dry
        self.m_tcs_percent = self.m_tcs/self.m_dry
        self.m_ttc_percent = self.m_ttc/self.m_dry
        return {
            "aocs_percent": self.m_aocs_percent,
            "cdh_percent": self.m_cdh_percent,
            "capture_system_percent": self.m_capture_system_percent,
            "eps_percent": self.m_eps_percent,
            "propulsion_percent": self.m_propulsion_percent,
            "structures_percent": self.m_structures_percent,
            "tcs_percent": self.m_tcs_percent,
            "ttc_percent": self.m_ttc_percent
        }