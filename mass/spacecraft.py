class SpaceCraft:
    """
    SpaceCraft digital twin class
    """

    def __init__(self):
        self.propulsion = None
        self.capture_system = None
        self.AOCS = None
        self.EPS = None
        self.structures = None
        self.TCS = None
        self.TTC = None
        self.CDH = None

        return

    def mass(self):
        """
        Compute the mass of the spacecraft

        """

        return