# polarizabilityCalc

A set of functions for calculating atomic polarizabilities given the atom's spectrum.

The code reads a .xlsx file that contains a list of the atomic transition lines in the following order:

Wavelength | Aki | Elem | ElemCode | En_low | J_low | En_up | J_up | Aki correction | Label | Ref

The explanation of each column is as follows:

1. Wavelength: the wavelength at which the transition occurs;
2. Aki: the spontaneous-emission Einstein Coefficient (in s^(-1));
3. Elem: The atomic number of the atom considered;
4. ElemCode: The name of the atom considered, along with its ionization state (e.g. HeI, HeII...);
5. En_low: the energy corresponding to the state from which the transition occurs;
6. J_low: the total angular momentum of the atom in the state from which the transition occurs;
7. En_up: the energy corresponding to the state up to which the transition occurs;
8. J_up: the total angular momentum of the atom in the state up to which the transition occurs;
9. Aki correction: 1 or 0 -- if the Aki was corrected for some reason. 
   For example: if the lines were theoretically calculated, and some of them had the Aki corrected for experimental values.
10. Label: if it's a "special" transition. For example: 583 nm intercombination line of ErI;
11. Ref: the reference from which the data on the transitions were taken.
