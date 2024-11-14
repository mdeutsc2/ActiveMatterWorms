module BoundaryTypes {
    enum BD_TYPE { CIRCLE, CARDIOID, EPICYCLOID, EPITROCHOID}

    record Boundary{
        var t: BD_TYPE;

        proc init(t: BD_TYPE) {
            this.t = t;
        }
        // proc init(t: string) {
            
        // }
    }

    proc boundary_type_init(input:int) {
        if (input == 1) {
            return new Boundary(t = BD_TYPE.CIRCLE);
        } else if (input == 2) {
            return new Boundary(t = BD_TYPE.CARDIOID);
        } else if (input == 3) {
            return new Boundary(t = BD_TYPE.EPICYCLOID);
        } else if (input == 4) {
            return new Boundary(t = BD_TYPE.EPITROCHOID);
        } else {
            writeln("invalid boundary type!");
            halt();
        }
    }
}