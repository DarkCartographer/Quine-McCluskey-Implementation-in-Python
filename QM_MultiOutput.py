import sys
import itertools
from collections import defaultdict


#Parse BLIF file to extract information
def parse_blif(file_content):
    lines = [line.strip() for line in file_content if line.strip() and not line.strip().startswith('#')]   #removes white-space and comments
    model_name = ""  #store model name
    inputs = []     #input variables
    outputs = []    #output variables
    names_blocks = []   #logic blocks
    i = 0

    #iterate through every line collecting information and storing them in the corresponding list
    while i < len(lines):
        line = lines[i]
        if line.startswith('.model'):
            model_name = line.split()[1]
            i += 1
        elif line.startswith('.inputs'):
            inputs.extend(line.split()[1:])
            i += 1
        elif line.startswith('.outputs'):
            outputs.extend(line.split()[1:])
            i += 1
        elif line.startswith('.names'):
            block_inputs_outputs = line.split()[1:]
            block = {'inputs_outputs': block_inputs_outputs, 'lines': []}
            i += 1
            while i < len(lines) and not lines[i].startswith('.'):
                block['lines'].append(lines[i])
                i += 1
            names_blocks.append(block)
        elif line.startswith('.end'):
            break
        else:
            i += 1
    return model_name, inputs, outputs, names_blocks

def simulate_blif(model_name, inputs, outputs, names_blocks):

    #track all variables with these
    variables = set(inputs)
    variables.update(outputs)   #add outputs
    intermediate_vars = set()   #intermideate variables (not a primary input or output)
    var_order = []  #keep track of the order to evaluate variables
    functions = {}


    #go through each logic block and build the corresponding truth table
    for block in names_blocks:
        block_vars = block['inputs_outputs']
        output_var = block_vars[-1]         #last variable is the output
        input_vars = block_vars[:-1]        #all others are inputs
        variables.update(block_vars)        #add them this set
        if output_var not in outputs:
            intermediate_vars.add(output_var)       #if not a primary output, mark it as intermediate
        var_order.append(output_var)
        # Build truth table for this block
        terms = []
        for line in block['lines']:
            if ' ' in line:
                in_part, out_part = line.split()
                if out_part == '1':
                    terms.append(in_part)
        functions[output_var] = {'inputs': input_vars, 'terms': terms}

    # Generate truth tables for all primary outputs
    input_vars = inputs
    n = len(input_vars)
    truth_tables = {output: {} for output in outputs}
    for bits in itertools.product([0, 1], repeat=n):
        assignment = dict(zip(input_vars, bits))

        # Evaluate intermediate variables
        for var in var_order:
            func = functions[var]
            input_vals = ''.join([str(assignment[v]) if v in assignment else '0' for v in func['inputs']])
            value = 0
            for term in func['terms']:
                match = True
                for idx, char in enumerate(term):
                    if char != '-' and char != input_vals[idx]:
                        match = False
                        break
                if match:
                    value = 1
                    break
            assignment[var] = value         #if a match was found, store the calculated value

        # Record the truth table for primary outputs
        input_key = ''.join(map(str, bits))     #convert the values to a binary string
        for output in outputs:
            output_value = assignment[output]
            truth_tables[output][input_key] = output_value
    return truth_tables, input_vars


#actually implement the QM method
def quine_mccluskey(truth_table, num_vars):
    minterms = [int(k, 2) for k, v in truth_table.items() if v == 1]
    minterms_bin = [format(m, f'0{num_vars}b') for m in minterms]

    # Group minterms by number of ones
    groups = defaultdict(list)
    for m in minterms_bin:
        ones = m.count('1')
        groups[ones].append(m)

    # Begin combining minterms
    prime_implicants = set()
    while True:
        new_groups = defaultdict(list)      # Dictionary for new groups after combination
        marked = set()
        keys = sorted(groups.keys())        #sort groups by number of ones they contain
        for i in range(len(keys)-1):        #iterate through adjacent groups
            for m1 in groups[keys[i]]:      #look through minterms in first group
                for m2 in groups[keys[i+1]]:#do the same for second group
                    diff = sum(c1 != c2 for c1, c2 in zip(m1, m2))  #count the differences in ones between the two terms
                    if diff == 1:       #if there is a difference of 1
                        #find the different position and ...
                        idx = next(j for j in range(len(m1)) if m1[j] != m2[j]) #create a new combined term
                        combined = m1[:idx] + '-' + m1[idx+1:]
                        new_groups[keys[i]].append(combined) #mark the original terms as combined
                        marked.add(m1)
                        marked.add(m2)

        #find unmarked terms and record them as prime implicants
        unmarked = set(itertools.chain.from_iterable(groups.values())) - marked
        prime_implicants.update(unmarked)
        if not new_groups:
            break
        groups = new_groups

    # Essential prime implicants
    chart = defaultdict(set)        #to keep track of covered minterms
    for m in minterms_bin:
        for p in prime_implicants:
            if match_minterm(p, m):     #check if the minterm is covered by the implicant
                chart[m].add(p)         #if it does, add it to the list for this term
    essential_prime_implicants = set()
    while chart:
        # Find primes that cover minterms uniquely
        essential = None
        for m, ps in chart.items():
            if len(ps) == 1:
                essential = next(iter(ps))
                break
        if essential:
            essential_prime_implicants.add(essential)
            # remove all minterms covered by this prime implicant
            to_remove = []
            for m in chart:
                if essential in chart[m]:
                    to_remove.append(m)
            for m in to_remove:
                del chart[m]
        else:
            # If no essential prime implicants, pick one with most minterms
            counts = defaultdict(int)
            for m, ps in chart.items():
                for p in ps:
                    counts[p] += 1
            if not counts:
                break
            essential = max(counts, key=counts.get)
            essential_prime_implicants.add(essential)
            # Remove all minterms covered by this prime implicant
            to_remove = []
            for m in chart:
                if essential in chart[m]:
                    to_remove.append(m)
            for m in to_remove:
                del chart[m]
    return essential_prime_implicants

def match_minterm(implicant, minterm):
    return all(ic == mc or ic == '-' for ic, mc in zip(implicant, minterm))

def main():
    # Prompt the user to enter the file name
    filename = input("Please enter the BLIF file name. Example <input.blif>: ").strip()

    try:
        # Try to open the file
        with open(filename, 'r') as f:
            file_content = f.readlines()
    except FileNotFoundError:
        # Handle the case where the file is not found
        print(f"Error: File '{filename}' not found. Please check the file name and try again.")
        return
    except Exception as e:
        # Handle other errors (e.g., permission issues)
        print(f"An error occurred while reading the file: {e}")
        return

    model_name, inputs, outputs, names_blocks = parse_blif(file_content)
    truth_tables, input_vars = simulate_blif(model_name, inputs, outputs, names_blocks)
    num_vars = len(input_vars)

    pla_lines = []
    pla_lines.append('.i {}'.format(num_vars))
    pla_lines.append('.o {}'.format(len(outputs)))
    pla_lines.append('.ilb {}'.format(' '.join(input_vars)))
    pla_lines.append('.ob {}'.format(' '.join(outputs)))

    # Dictionary to store PLA terms: key is implicant, value is output pattern
    pla_terms_dict = {}

    for idx, output in enumerate(outputs):
        output_truth_table = {k: v for k, v in truth_tables[output].items()}
        essential_implicants = quine_mccluskey(output_truth_table, num_vars)
        for implicant in essential_implicants:
            if implicant in pla_terms_dict:
                # Update existing output pattern
                output_values = list(pla_terms_dict[implicant])
                output_values[idx] = '1'
                pla_terms_dict[implicant] = ''.join(output_values)
            else:
                # Create new output pattern
                output_values = ['0'] * len(outputs)
                output_values[idx] = '1'
                pla_terms_dict[implicant] = ''.join(output_values)

    # Prepare the PLA terms
    pla_terms = ['{} {}'.format(implicant, output_values) for implicant, output_values in pla_terms_dict.items()]
    pla_lines.append('.p {}'.format(len(pla_terms)))
    pla_lines.extend(pla_terms)
    pla_lines.append('.e')
    pla_output = '\n'.join(pla_lines)
    print(pla_output)

if __name__ == '__main__':
    main()
