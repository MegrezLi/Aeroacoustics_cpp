"""One-time, reviewable AST translation of the already validated Python kernels.
Produces ordinary C++ source; no Python or interpreter at library runtime.
"""
from pathlib import Path
import ast
ROOT=Path(__file__).resolve().parents[1]
tree=ast.parse((ROOT/'reference/python/aeroacoustics.py').read_text(encoding='utf-8'))
class ReservedNames(ast.NodeTransformer):
    def visit_Name(self,node):
        if node.id=='switch':node.id='separated'
        return node
tree=ReservedNames().visit(tree)
NAMES='log10aa lblvs tblte tipnois inflownoise blunt g5comp amin amax bmin bmax a0comp thick directh_te directh_le directl simple_guidati'.split()
ARRAY_FUNCS={'lblvs','tipnois','inflownoise','blunt','simple_guidati'}
MATH={'sqrt','log10','log','exp','sin','cos','tanh','abs'}

def expr(n):
    if isinstance(n,ast.Name):return n.id
    if isinstance(n,ast.Constant):
        if isinstance(n.value,bool):return 'true' if n.value else 'false'
        return repr(n.value)
    if isinstance(n,ast.Attribute):return expr(n.value)+'.'+n.attr
    if isinstance(n,ast.BinOp):
        a,b=expr(n.left),expr(n.right)
        if isinstance(n.op,ast.Pow):return f'std::pow({a}, {b})'
        op={ast.Add:'+',ast.Sub:'-',ast.Mult:'*',ast.Div:'/',ast.Mod:'%'}[type(n.op)]
        # Force floating division for Python literal integer ratios (e.g. 7/3).
        if isinstance(n.op,ast.Div):return f'({a} / static_cast<double>({b}))'
        return f'({a} {op} {b})'
    if isinstance(n,ast.UnaryOp):return {ast.USub:'-',ast.UAdd:'+',ast.Not:'!'}[type(n.op)]+'('+expr(n.operand)+')'
    if isinstance(n,ast.BoolOp):return '('+(' && ' if isinstance(n.op,ast.And) else ' || ').join(map(expr,n.values))+')'
    if isinstance(n,ast.Compare):
        ops={ast.Eq:'==',ast.NotEq:'!=',ast.Lt:'<',ast.LtE:'<=',ast.Gt:'>',ast.GtE:'>='}
        return '('+expr(n.left)+' '+ops[type(n.ops[0])]+' '+expr(n.comparators[0])+')'
    if isinstance(n,ast.Call):
        fn=expr(n.func);args=list(map(expr,n.args))
        if fn=='len':return f'{args[0]}.size()'
        if fn in ('max','min'):fn='std::'+fn+'<double>'
        if fn.startswith('np.'):fn='std::'+fn[3:]
        if fn=='abs':fn='std::abs'
        return fn+'('+', '.join(args)+')'
    if isinstance(n,ast.Subscript):return expr(n.value)+'['+expr(n.slice)+']'
    if isinstance(n,ast.Tuple):return '{'+', '.join(map(expr,n.elts))+'}'
    raise ValueError(ast.dump(n))

defs=[];decls=[]
for f in tree.body:
    if not isinstance(f,ast.FunctionDef) or f.name not in NAMES:continue
    args=[x.arg for x in f.args.args]
    ret='Spectrum' if f.name in ARRAY_FUNCS else 'double'
    if f.name=='tblte':ret='std::tuple<Spectrum, Spectrum, Spectrum>'
    if f.name=='thick':ret='std::tuple<double, double, double>'
    sig=ret+' '+f.name+'('+', '.join(('const Parameters& p' if a=='p' else 'double '+a) for a in args)+')'
    decls.append(sig+';')
    arrays={};variables=set()
    for n in ast.walk(f):
        if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Store):variables.add(n.id)
        if isinstance(n,ast.Assign) and isinstance(n.value,ast.Call) and expr(n.value.func)=='np.zeros':
            arrays[n.targets[0].id]=expr(n.value.args[0])
    body=[sig+' {']
    loopvars={n.target.id for n in ast.walk(f) if isinstance(n,ast.For)}
    for v in sorted(variables-set(args)-set(arrays)):
        body.append('    '+('int' if v in loopvars else 'double')+' '+v+' = 0;')
    for a,size in arrays.items():body.append(f'    Spectrum {a}({size}, 0.);')
    def emit(n,depth=1):
        tab='    '*depth
        if isinstance(n,ast.Expr) and isinstance(n.value,ast.Constant):return
        if isinstance(n,ast.Pass):body.append(tab+';');return
        if isinstance(n,ast.Assign):
            if isinstance(n.value,ast.Call) and expr(n.value.func)=='np.zeros':return
            target=n.targets[0]
            if isinstance(target,ast.Tuple):
                lhs='std::tie('+', '.join(map(expr,target.elts))+')'
            elif isinstance(target,ast.Subscript) and isinstance(target.slice,ast.Slice):
                v=expr(target.value);body.append(f'{tab}std::fill({v}.begin(), {v}.end(), {expr(n.value)});');return
            else:lhs=expr(target)
            body.append(tab+lhs+' = '+expr(n.value)+';');return
        if isinstance(n,ast.If):
            body.append(tab+'if ('+expr(n.test)+') {')
            for sub in n.body:emit(sub,depth+1)
            if n.orelse:
                body.append(tab+'} else {')
                for sub in n.orelse:emit(sub,depth+1)
            body.append(tab+'}');return
        if isinstance(n,ast.For):
            bounds=n.iter.args
            var=expr(n.target)
            body.append(tab+f'for ({var} = {expr(bounds[0])}; {var} < {expr(bounds[1])}; ++{var}) {{')
            for sub in n.body:emit(sub,depth+1)
            body.append(tab+'}');return
        if isinstance(n,ast.Return):body.append(tab+'return '+expr(n.value)+';');return
        raise ValueError(ast.dump(n))
    for n in f.body:emit(n)
    body.append('}')
    defs.append('\n'.join(body))
(ROOT/'include/kernels.hpp').write_text('#pragma once\n#include "aeroacoustics.hpp"\nnamespace aeroacoustics {\n'+'\n'.join(decls)+'\n}\n')
source='''// Numerical expressions ported from OpenFAST (Apache-2.0), see NOTICE.
// Materialized by tools/generate_kernels.py; no interpreter at runtime.
#include "kernels.hpp"
#include <algorithm>
#include <cmath>
namespace aeroacoustics {
constexpr double aa_epsilon=1e-16, pi=3.14159265358979323846, twopi=2*pi;
constexpr int x_blmethod_tables=2, itrip_none=0, itrip_heavy=1, itrip_light=2;
'''
(ROOT/'src/kernels.cpp').write_text(source+'\n\n'.join(defs)+'\n}\n')
# Copy the exact 61-point Kronrod weights, without inventing a new quadrature.
with (ROOT/'include/quadrature_data.hpp').open('w') as out:
    out.write('#pragma once\n#include <array>\nnamespace aeroacoustics {\n')
    for n in tree.body:
        if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ('_XGK','_WGK','_WG'):
            values=ast.literal_eval(n.value.args[0]);name=n.targets[0].id[1:].lower()
            out.write(f'inline constexpr std::array<double,{len(values)}> {name} = '+ '{'+', '.join(map(repr,values))+'};\n')
    out.write('}\n')
print('Generated',len(defs),'C++ kernel functions')
