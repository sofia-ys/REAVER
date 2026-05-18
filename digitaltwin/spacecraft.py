from subsystems import Propulsion, AOCS, Structures, TCS, CaptureSystem, EPS, TTC, CDH

class SpaceCraft:
    """
    SpaceCraft digital twin class
    """

    def __init__(self, n_targets, dv_list, m_debris_list, Isp, sc_type, n_redundancy, capture_type, mission_type):  # parameters needed for all the subsystem classes here
        self.n_targets = n_targets
        self.dv_list = dv_list
        self.m_debris_list = m_debris_list
        self.Isp = Isp
        self.sc_type = sc_type
        self.n_redundancy = n_redundancy
        self.capture_type = capture_type
        self.mission_type = mission_type

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

        return m_dry