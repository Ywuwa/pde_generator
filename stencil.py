# stencil.py
from collections import defaultdict
from ast_nodes import Const, Add, Mul
from symbolic import simplify_ast

class StencilExpr:
    def __init__(self):
        # (dx,dy,dz) -> coef (AST)
        self.terms = defaultdict(lambda: Const(0))

    def copy(self):
        s = StencilExpr()
        for k, v in self.terms.items():
            s.terms[k] = v
        return s

    def add_term(self, offset, coef):
        """
        offset: (dx,dy,dz)
        coef: AST (Const, Mul, ...)
        """
        if coef is None:
            return

        if offset in self.terms:
            self.terms[offset] = Add(self.terms[offset], coef)
        else:
            self.terms[offset] = coef

    def __add__(self, other):
        result = self.copy()
        for offset, coef in other.terms.items():
            result.add_term(offset, coef)
        return result.simplify()

    def scale(self, scalar):
        result = StencilExpr()
        for offset, coef in self.terms.items():
            result.terms[offset] = Mul(scalar, coef)
        return result.simplify()

    def __mul__(self, other):
        """
        Свёртка stencil (ключевая операция!)
        """
        result = StencilExpr()
        for (dx1,dy1,dz1), c1 in self.terms.items():
            for (dx2,dy2,dz2), c2 in other.terms.items():
                offset = (dx1+dx2, dy1+dy2, dz1+dz2)
                coef = Mul(c1, c2)
                result.add_term(offset, coef)
        return result.simplify()
      
    def __repr__(self):
      parts = []
      for (dx,dy,dz), coef in sorted(self.terms.items()):
          parts.append(f"{coef} @ ({dx},{dy},{dz})")
      return " + ".join(parts) if parts else "0"
    
    def simplify(self):
      """
      Упрощает AST-выражения (убирает нули, обнуляющиеся слагаемые и т.д.)
      """
      cleaned = {}
      for k, v in self.terms.items():
          v = simplify_ast(v)
          if isinstance(v, Const) and v.value == 0:
              continue
          cleaned[k] = v
  
      self.terms = cleaned
      return self