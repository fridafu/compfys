#ifndef Nproblem_h
#define Nproblem_h

template <class T>
class Nproblem
{
    private:
        T objects[];
        int numobj;
        vec t;
        double dt;
        vec accel1(numobj);
        double accel2;

    public:
        Nproblem(T);
        double solve(vec);
}

class Verlet : public Nproblem
{
    public:
        void advance();
}

#endif
