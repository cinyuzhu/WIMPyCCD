class DMParams:
    def __init__(self,
                 m_D_GeV = 5,
                 sigma_pb = 1,
                 rho_GeV_cm3 = 0.3,
                 v0_kms = 238,
                 vE_kms = 263,
                 vesc_kms= 544,
                 confidence_level = 0.9,
                 mass_number_A = 28,
                 ):
        
        # Natural units
        self.m_D_GeV = m_D_GeV
        self.sigma_pb = sigma_pb
        self.rho_GeV_cm3 = rho_GeV_cm3
        self.v0_kms = v0_kms
        self.vE_kms = vE_kms
        self.vesc_kms = vesc_kms
        self.confidence_level = confidence_level
        self.mass_number_A = mass_number_A
        self.N0 = 6.02e26
    
    # to keep the following parameters updated in real time
    @property
    def m_T_GeV(self):
        return self.mass_number_A * 0.938
        
    @property
    def m_D_kg(self):
        return self.m_D_GeV * 1.78e-27

    @property
    def m_T_kg(self):
        return self.m_T_GeV * 1.78e-27

    @property
    def mu_n(self):
        return self.m_D_kg * self.m_T_kg / (self.m_D_kg + self.m_T_kg)

    @property
    def sigma_m2(self):
        return self.sigma_pb * 1e-40

    @property
    def rho_SI(self):
        return self.rho_GeV_cm3 * 1.783e-21

    @property
    def v0_mps(self):
        return self.v0_kms * 1e3

    @property
    def vE_mps(self):
        return self.vE_kms * 1e3

    @property
    def vesc_mps(self):
        return self.vesc_kms * 1e3
