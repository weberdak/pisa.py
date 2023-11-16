import numpy


class Wheel:
    def __init__(self, tau=20.0, rho=0.0, order=1.0, phi=-63.0, psi=-42.0,
                 beta=17.0, dc=10.735, pas=(57.3, 81.2, 228.1), seq='',
                 seq_start=1, rho_start=1, nres=18, flip=0.0, period=3.6,
                 negative_dc=True, aCaCN=117.5, aCNCa=124.0, aNCaC=107.4,
                 aCaNH=116.0, bCaC=1.52, bCN=1.35, bNCa=1.45):
        self.tau = 20.0
        self.rho = 0.0
        self.order = 1.0
        self.phi = -63.0
        self.psi = -42.0
        self.beta = 17.0
        self.dc = 10.735
        self.pas = (57.3, 81.2, 228.1)
        self.seq = ''
        self.seq_start = 1
        self.rho_start = 1
        self.nres = 18
        self.flip = 0.0
        self.period = 3.6
        self.negative_dc = True
        self.aCaCN = 117.5
        self.aCNCa = 124.0
        self.aNCaC = 107.4
        self.aCaNH = 116.0
        self.bCaC = 1.52
        self.bCN = 1.35
        self.bNCa = 1.45

    @property
    def tau(self):
        return numpy.degrees(self._tau)
    
    @tau.setter
    def tau(self, value):
        if value < 0 or value > 90:
            raise ValueError('Tau must be between 0 and 90 degrees')
        self._tau = numpy.radians(value)

    @property
    def rho(self):
        return numpy.degrees(self._rho)
    
    @rho.setter
    def rho(self, value):
        if value < 0 or value > 360:
            raise ValueError('Rho must be between 0 and 360 degrees')
        self._rho = numpy.radians(value)
    
    @property
    def order(self):
        return self._order
    
    @order.setter
    def order(self, value):
        if value < 0 or value > 1:
            raise ValueError('Order must be between 0 and 1')
        self._order = value
    
    @property
    def phi(self):
        return numpy.degrees(self._phi)
    
    @phi.setter
    def phi(self, value):
        if value < -180 or value > 180:
            raise ValueError('Phi must be between -180 and 180 degrees')
        self._phi = numpy.radians(value)
    
    @property
    def psi(self):
        return numpy.degrees(self._psi)
    
    @psi.setter
    def psi(self, value):
        if value < -180 or value > 180:
            raise ValueError('Psi must be between -180 and 180 degrees')
        self._psi = numpy.radians(value)
    
    @property
    def beta(self):
        return numpy.degrees(self._beta)
    
    @beta.setter
    def beta(self, value):
        if value < -180 or value > 180:
            raise ValueError('Beta must be between -180 and 180 degrees')
        self._beta = numpy.radians(value)
        
    @property
    def dc(self):
        return self._dc
    
    @dc.setter
    def dc(self, value):
        if value < 0:
            raise ValueError('DC must be greater than 0')
        self._dc = value
        
    @property
    def pas(self):
        return self._pas
    
    @pas.setter
    def pas(self, value):
        if len(value) != 3:
            raise ValueError('PAS must be a tuple of length 3')
        self._pas = value
        
    @property
    def seq(self):
        return self._seq
    
    @seq.setter
    def seq(self, value):
        if not value.isupper():
            raise ValueError('Sequence must be uppercase letters')
        self._nres = len(value)
        self._seq = value
        
    @property
    def seq_start(self):
        return self._seq_start
    
    @seq_start.setter
    def seq_start(self, value):
        if value < 1:
            raise ValueError('Sequence start must be greater than 0')
        self._seq_start = value