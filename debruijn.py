# from Bio import SeqIO
# import argparse
# import matplotlib.pyplot as plt
# from matplotlib.patches import Circle, FancyArrowPatch
# import math
# from matplotlib.backends.backend_pdf import PdfPages
# from collections import defaultdict, Counter

# def build_debruijn_graph(kmers):
#     """
#     Build a de Bruijn graph from k-mers and their counts
#     Returns edges dictionary and nodes set
#     """
#     edges = defaultdict(Counter)
#     nodes = set()
    
#     for kmer in kmers:
#         prefix = kmer[:-1]
#         suffix = kmer[1:]
#         nodes.add(prefix)
#         nodes.add(suffix)
#         edges[prefix][suffix] += 1
        
#     return edges, nodes

# def find_eulerian_path(edges, nodes):
#     """Find Eulerian path preserving multiplicities"""
#     if not nodes:
#         return []
        
#     # Find start node (any node with outgoing edges)
#     start = next(iter(edges.keys()))
    
#     path = []
#     stack = [start]
#     edges_remaining = {node: Counter(edges[node]) for node in edges}
    
#     while stack:
#         current = stack[-1]
        
#         # If current node has remaining outgoing edges
#         if edges_remaining[current]:
#             # Get next node that still has remaining edges
#             next_node = next(iter(edges_remaining[current]))
            
#             # Add k-mer to path and decrement edge count
#             path.append(current + next_node[-1])
#             edges_remaining[current][next_node] -= 1
#             if edges_remaining[current][next_node] == 0:
#                 del edges_remaining[current][next_node]
            
#             stack.append(next_node)
#         else:
#             stack.pop()
            
#     return path

# def find_kmers(sequence, k):
#     """
#     Find k-mers in sequence order
#     Returns list of k-mers in order of appearance
#     """
#     # Add k-1 characters to handle circular nature
#     extended_sequence = sequence + sequence[:k-1]
    
#     # Get k-mers in order
#     kmers = []
#     for i in range(len(sequence)):
#         kmer = extended_sequence[i:i+k]
#         kmers.append(kmer)
    
#     return kmers

# def get_nodes(kmers):
#     """
#     Find the nodes, which are (k-1)-mers
#     Returns list of unique (k-1)-mers
#     """
#     nodes = set()
#     for kmer in kmers:
#         nodes.add(kmer[:-1])
#         nodes.add(kmer[1:])
#     return sorted(list(nodes))

# def get_node_positions(nodes):
#     """
#     Calculate node positions in a circle
#     Returns dictionary mapping nodes to (x,y) coordinates
#     """
#     positions = {}
#     num_nodes = len(nodes)
#     radius = 1  # Unit circle
    
#     for i in range(num_nodes):
#         # Calculate position on circle
#         angle = 2 * math.pi * i / num_nodes
#         x = radius * math.cos(angle)
#         y = radius * math.sin(angle)
#         positions[nodes[i]] = (x, y)
    
#     return positions

# def create_figure():
#     """Create and return figure and axis"""
#     fig, ax = plt.subplots(figsize=(12, 12))
#     ax.set_aspect('equal')
#     ax.set_xlim(-1.5, 1.5)
#     ax.set_ylim(-1.5, 1.5)
#     plt.axis('off')
#     return fig, ax

# def draw_nodes(ax, nodes, positions):
#     """Draw nodes as circles with labels"""
#     node_radius = 0.1
#     for node, (x, y) in positions.items():
#         circle = Circle((x, y), node_radius, 
#                       facecolor='lightblue',
#                       edgecolor='black',
#                       alpha=0.6)
#         ax.add_patch(circle)
        
#         # Add node label
#         plt.text(x, y, node,
#                 horizontalalignment='center',
#                 verticalalignment='center',
#                 fontweight='bold')
#     return node_radius

# def draw_edges(ax, kmers, positions, node_radius):
#     """Draw edges as arrows with labels showing true sequence order"""
#     # Create a mapping of edge to count and track order
#     edge_counts = defaultdict(int)
#     edge_orders = []  # List to maintain sequence order of edges
    
#     # First pass - count edges and record order
#     for kmer in kmers:
#         edge = (kmer[:-1], kmer[1:])
#         edge_counts[edge] += 1
#         edge_orders.append(edge)
    
#     # Track which edges we've drawn
#     drawn_edges = set()
    
#     # Draw each unique edge once
#     edge_number = 1
#     for edge in edge_counts:
#         if edge in drawn_edges:
#             continue
            
#         start_node, end_node = edge
#         count = edge_counts[edge]
#         start_pos = positions[start_node]
#         end_pos = positions[end_node]
        
#         # Calculate actual start and end points
#         dx = end_pos[0] - start_pos[0]
#         dy = end_pos[1] - start_pos[1]
#         length = math.sqrt(dx*dx + dy*dy)
        
#         # Skip if length is 0 (same position nodes)
#         if length == 0:
#             continue

#         # Adjust start and end points
#         start_x = start_pos[0] + (dx * node_radius / length)
#         start_y = start_pos[1] + (dy * node_radius / length)
#         end_x = end_pos[0] - (dx * node_radius / length)
#         end_y = end_pos[1] - (dy * node_radius / length)
        
#         # Draw arrow
#         arrow = FancyArrowPatch(
#             (start_x, start_y), (end_x, end_y),
#             arrowstyle='-|>',
#             mutation_scale=20,
#             linewidth=1,
#             color='gray'
#         )
#         ax.add_patch(arrow)
        
#         # Find all positions where this edge occurs
#         edge_positions = []
#         for i, e in enumerate(edge_orders):
#             if e == edge:
#                 edge_positions.append(i + 1)  # +1 for 1-based numbering
        
#         # Draw labels
#         kmer = start_node + end_node[-1]
#         for i, pos in enumerate(edge_positions):
#             # Calculate position for this label
#             frac = (i + 1)/(count + 1)  # Spread labels evenly
#             mid_x = start_x + (end_x - start_x) * frac
#             mid_y = start_y + (end_y - start_y) * frac
            
#             if count > 1:
#                 label = f"{pos}:{kmer}(x{count}, #{i+1})"
#             else:
#                 label = f"{pos}:{kmer}"
                
#             plt.text(mid_x, mid_y, label,
#                     bbox=dict(facecolor='white', edgecolor='black', alpha=0.7, pad=2),
#                     horizontalalignment='center',
#                     verticalalignment='center')
        
#         drawn_edges.add(edge)

# def assemble_sequence(kmers, original_length):
#     """
#     Assemble sequence from k-mers in their original order
#     Returns assembled sequence of specified length
#     """
#     if not kmers:
#         return None
        
#     # For first k-mer, take all but last character
#     sequence = kmers[0][:-1]
    
#     # For each k-mer, take last character
#     for kmer in kmers:
#         sequence += kmer[-1]
        
#     return sequence[:original_length]

# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('fasta_file')
#     parser.add_argument('-k', type=int)
#     args = parser.parse_args()
#     k = args.k
    
#     # Read sequence
#     sequence = str(next(SeqIO.parse(args.fasta_file, "fasta")).seq)
    
#     if len(sequence) < k:
#         print(f"Error: Sequence length ({len(sequence)}) must be >= k ({k})")
#         return
    
#     # Find kmers in sequence order
#     kmers = find_kmers(sequence, k)
#     kmer_counts = Counter(kmers)
    
#     print("\nk-mers (edges) with counts:")
#     for kmer in kmers:  # Print in sequence order
#         print(f"{kmer}: appears {kmer_counts[kmer]}x")
    
#     # Get nodes (k-1)mers
#     nodes = get_nodes(kmers)
#     print("\n(k-1)-mers (nodes):")
#     for i, node in enumerate(nodes, 1):
#         print(f"{i}. {node}")
    
#     print("\nEdges in sequence order:")
#     for i, kmer in enumerate(kmers, 1):
#         print(f"{i}. {kmer[:-1]} -> {kmer[1:]} ({kmer})")
    
#     # Calculate node positions
#     positions = get_node_positions(nodes)
    
#     # Create figure
#     fig, ax = create_figure()
    
#     # Draw nodes
#     node_radius = draw_nodes(ax, nodes, positions)
    
#     # Draw edges with multiplicities
#     draw_edges(ax, kmers, positions, node_radius)
    
#     # Create pdf file for output
#     with PdfPages('debruijn_visual.pdf') as pdf:
#         pdf.savefig(fig, bbox_inches='tight')
#     plt.close()
    
#     # Assemble sequence
#     assembled = assemble_sequence(kmers, len(sequence))
    
#     # Print results
#     print("\nDe Bruijn Graph Assembly Results:")
#     print(f"k-mer size: {k}")
#     print(f"Number of unique k-mers: {len(kmer_counts)}")
#     print(f"Total k-mer count (with multiplicities): {len(kmers)}")
#     print(f"Number of nodes: {len(nodes)}")
#     print(f"Original sequence: {sequence}")
#     print(f"Assembled sequence: {assembled}")
#     print(f"Assembly matches original: {sequence == assembled}")

#     with open('output.txt', 'w') as outfile:
#         outfile.write(f"De Bruijn Graph Assembly Results:\n")
#         outfile.write(f"k-mer size: {k}\n")
#         outfile.write(f"Number of unique k-mers: {len(kmer_counts)}\n")
#         outfile.write(f"Total k-mer count (with multiplicities): {len(kmers)}\n")
#         outfile.write(f"Number of nodes: {len(nodes)}\n")
#         outfile.write(f"Original sequence: {sequence}\n")
#         outfile.write(f"Assembled sequence: {assembled}\n")
#         outfile.write(f"Assembly matches original: {sequence == assembled}\n")

# if __name__ == "__main__":
#     main()



# from Bio import SeqIO
# import argparse
# import matplotlib.pyplot as plt
# from matplotlib.patches import Circle, FancyArrowPatch
# import math
# from matplotlib.backends.backend_pdf import PdfPages
# from collections import defaultdict, Counter
# import random

# def introduce_errors(sequence, error_rate=0.05, min_gap_size=3):
#     """
#     Introduce errors in the sequence to create gaps/contigs
#     """
#     seq_list = list(sequence)
#     error_positions = []
    
#     # Introduce gaps (removing nucleotides)
#     i = 0
#     while i < len(seq_list):
#         if random.random() < error_rate:
#             gap_size = random.randint(min_gap_size, min_gap_size + 2)
#             if i + gap_size < len(seq_list):
#                 for j in range(gap_size):
#                     seq_list[i + j] = 'N'
#                 error_positions.extend(range(i, i + gap_size))
#                 i += gap_size
#             else:
#                 break
#         i += 1
    
#     return ''.join(seq_list), error_positions

# def find_strongly_connected_components(edges, nodes):
#     """
#     Find strongly connected components in the graph using Kosaraju's algorithm
#     """
#     def dfs_forward(node, visited, order):
#         visited.add(node)
#         if node in edges:
#             for next_node in edges[node]:
#                 if next_node not in visited:
#                     dfs_forward(next_node, visited, order)
#         order.append(node)

#     def dfs_reverse(node, visited, component):
#         visited.add(node)
#         component.append(node)
#         for prev_node in nodes:
#             if prev_node in edges and node in edges[prev_node] and prev_node not in visited:
#                 dfs_reverse(prev_node, visited, component)

#     # First DFS to get order
#     visited = set()
#     order = []
#     for node in nodes:
#         if node not in visited:
#             dfs_forward(node, visited, order)

#     # Second DFS to find SCCs
#     visited = set()
#     components = []
#     for node in reversed(order):
#         if node not in visited:
#             component = []
#             dfs_reverse(node, visited, component)
#             if component:
#                 components.append(component)

#     return components

# def find_contigs(edges, nodes):
#     """
#     Find contigs in the De Bruijn graph by identifying maximal paths
#     """
#     components = find_strongly_connected_components(edges, nodes)
#     contigs = []
    
#     for component in components:
#         # For each component, find maximal non-branching paths
#         visited_edges = set()
        
#         for start_node in component:
#             if start_node not in edges:
#                 continue
                
#             # Check if this node is the start of a non-branching path
#             while start_node in edges and len(edges[start_node]) == 1:
#                 current_path = [start_node]
#                 current = start_node
                
#                 # Follow the path until we reach a branch or visited edge
#                 while current in edges and len(edges[current]) == 1:
#                     next_node = next(iter(edges[current]))
#                     edge = (current, next_node)
                    
#                     if edge in visited_edges:
#                         break
                        
#                     visited_edges.add(edge)
#                     current_path.append(next_node)
#                     current = next_node
                    
#                     # Break if we've returned to start or reached a branch
#                     if current == start_node or (current in edges and len(edges[current]) != 1):
#                         break
                
#                 if len(current_path) > 1:
#                     contigs.append(current_path)
#                 break
    
#     return contigs

# def assemble_contigs(contigs, kmers):
#     """
#     Assemble sequences for each contig
#     """
#     assembled_contigs = []
    
#     for contig in contigs:
#         if len(contig) < 2:  # Skip single-node contigs
#             continue
            
#         # Start with first node
#         sequence = contig[0]
        
#         # Add last character of each subsequent node
#         for node in contig[1:]:
#             sequence += node[-1]
            
#         assembled_contigs.append(sequence)
    
#     return assembled_contigs

# def get_longest_contig(assembled_contigs):
#     """
#     Find the longest contig from the assembled sequences
    
#     Parameters:
#     assembled_contigs (list): List of assembled contig sequences
    
#     Returns:
#     tuple: (longest_contig, length) or (None, 0) if no contigs
#     """
#     if not assembled_contigs:
#         return None, 0
        
#     longest_contig = max(assembled_contigs, key=len)
#     return longest_contig, len(longest_contig)

# def build_debruijn_graph(kmers):
#     """
#     Build a de Bruijn graph from k-mers and their counts
#     """
#     edges = defaultdict(Counter)
#     nodes = set()
    
#     for kmer in kmers:
#         # Skip k-mers containing errors (N's)
#         if 'N' in kmer:
#             continue
            
#         prefix = kmer[:-1]
#         suffix = kmer[1:]
#         nodes.add(prefix)
#         nodes.add(suffix)
#         edges[prefix][suffix] += 1
        
#     return edges, nodes

# def find_kmers(sequence, k):
#     """
#     Find k-mers in sequence order, skipping k-mers containing errors
#     """
#     # Add k-1 characters to handle circular nature
#     extended_sequence = sequence + sequence[:k-1]
    
#     # Get k-mers in order
#     kmers = []
#     for i in range(len(sequence)):
#         kmer = extended_sequence[i:i+k]
#         # Only include k-mers that don't contain errors
#         if 'N' not in kmer:
#             kmers.append(kmer)
    
#     return kmers

# def get_node_positions(nodes):
#     """Calculate node positions in a circle"""
#     positions = {}
#     num_nodes = len(nodes)
#     radius = 1  # Unit circle
    
#     for i, node in enumerate(nodes):
#         angle = 2 * math.pi * i / num_nodes
#         x = radius * math.cos(angle)
#         y = radius * math.sin(angle)
#         positions[node] = (x, y)
    
#     return positions

# def create_figure():
#     """Create and return figure and axis"""
#     fig, ax = plt.subplots(figsize=(12, 12))
#     ax.set_aspect('equal')
#     ax.set_xlim(-1.5, 1.5)
#     ax.set_ylim(-1.5, 1.5)
#     plt.axis('off')
#     return fig, ax

# def draw_nodes(ax, nodes, positions):
#     """Draw nodes as circles with labels"""
#     node_radius = 0.1
#     for node, (x, y) in positions.items():
#         circle = Circle((x, y), node_radius, 
#                       facecolor='lightblue',
#                       edgecolor='black',
#                       alpha=0.6)
#         ax.add_patch(circle)
#         plt.text(x, y, node,
#                 horizontalalignment='center',
#                 verticalalignment='center',
#                 fontweight='bold')
#     return node_radius

# def draw_edges(ax, kmers, positions, node_radius):
#     """Draw edges as arrows with labels"""
#     edge_counts = defaultdict(int)
#     edge_orders = []
    
#     for kmer in kmers:
#         if 'N' not in kmer:  # Only process valid k-mers
#             edge = (kmer[:-1], kmer[1:])
#             edge_counts[edge] += 1
#             edge_orders.append(edge)
    
#     drawn_edges = set()
    
#     for edge in edge_counts:
#         if edge in drawn_edges:
#             continue
            
#         start_node, end_node = edge
#         count = edge_counts[edge]
        
#         if start_node not in positions or end_node not in positions:
#             continue
            
#         start_pos = positions[start_node]
#         end_pos = positions[end_node]
        
#         dx = end_pos[0] - start_pos[0]
#         dy = end_pos[1] - start_pos[1]
#         length = math.sqrt(dx*dx + dy*dy)
        
#         if length == 0:
#             continue

#         start_x = start_pos[0] + (dx * node_radius / length)
#         start_y = start_pos[1] + (dy * node_radius / length)
#         end_x = end_pos[0] - (dx * node_radius / length)
#         end_y = end_pos[1] - (dy * node_radius / length)
        
#         arrow = FancyArrowPatch(
#             (start_x, start_y), (end_x, end_y),
#             arrowstyle='-|>',
#             mutation_scale=20,
#             linewidth=1,
#             color='gray'
#         )
#         ax.add_patch(arrow)
        
#         edge_positions = []
#         for i, e in enumerate(edge_orders):
#             if e == edge:
#                 edge_positions.append(i + 1)
        
#         kmer = start_node + end_node[-1]
#         for i, pos in enumerate(edge_positions):
#             frac = (i + 1)/(count + 1)
#             mid_x = start_x + (end_x - start_x) * frac
#             mid_y = start_y + (end_y - start_y) * frac
            
#             label = f"{pos}:{kmer}"
#             if count > 1:
#                 label += f"(x{count}, #{i+1})"
                
#             plt.text(mid_x, mid_y, label,
#                     bbox=dict(facecolor='white', edgecolor='black', alpha=0.7, pad=2),
#                     horizontalalignment='center',
#                     verticalalignment='center')
        
#         drawn_edges.add(edge)

# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('fasta_file')
#     parser.add_argument('-k', type=int)
#     parser.add_argument('--error_rate', type=float, default=0.05,
#                       help='Rate of error introduction (default: 0.05)')
#     args = parser.parse_args()
#     k = args.k
    
#     # Read sequence
#     sequence = str(next(SeqIO.parse(args.fasta_file, "fasta")).seq)
    
#     if len(sequence) < k:
#         print(f"Error: Sequence length ({len(sequence)}) must be >= k ({k})")
#         return
    
#     # Introduce errors
#     sequence_with_errors, error_positions = introduce_errors(sequence, args.error_rate)
    
#     # Find kmers in sequence order
#     kmers = find_kmers(sequence_with_errors, k)
#     kmer_counts = Counter(kmers)
    
#     # Build de Bruijn graph
#     edges, nodes = build_debruijn_graph(kmers)
    
#     # Find contigs
#     contigs = find_contigs(edges, nodes)
#     assembled_contigs = assemble_contigs(contigs, kmers)
    
#     print("\nDe Bruijn Graph Assembly Results with Errors:")
#     print(f"k-mer size: {k}")
#     print(f"Error rate: {args.error_rate}")
#     print(f"Number of error positions: {len(error_positions)}")
#     print(f"Number of contigs: {len(contigs)}")
#     print(f"\nOriginal sequence length: {len(sequence)}")
#     print(f"Number of k-mers: {len(kmers)}")
#     print(f"Number of nodes: {len(nodes)}")
    
#     print("\nContigs:")
#     for i, contig in enumerate(assembled_contigs, 1):
#         print(f"Contig {i} (length {len(contig)}):")
#         print(contig)
    
#     # Get and print longest contig
#     longest_contig, longest_length = get_longest_contig(assembled_contigs)
#     if longest_contig:
#         print(f"\nLongest contig (length {longest_length}):")
#         print(longest_contig)
#     else:
#         print("\nNo contigs found")
    
#     # Calculate node positions and create visualization
#     positions = get_node_positions(list(nodes))
#     fig, ax = create_figure()
#     node_radius = draw_nodes(ax, nodes, positions)
#     draw_edges(ax, kmers, positions, node_radius)
    
#     # Save visualization
#     with PdfPages('debruijn_visual_with_errors.pdf') as pdf:
#         pdf.savefig(fig, bbox_inches='tight')
#     plt.close()
    
#     # Save detailed results to file
#     with open('output_with_errors.txt', 'w') as outfile:
#         outfile.write(f"De Bruijn Graph Assembly Results with Errors:\n")
#         outfile.write(f"k-mer size: {k}\n")
#         outfile.write(f"Error rate: {args.error_rate}\n")
#         outfile.write(f"Number of error positions: {len(error_positions)}\n")
#         outfile.write(f"Error positions: {error_positions}\n")
#         outfile.write(f"Original sequence: {sequence}\n")
#         outfile.write(f"Sequence with errors: {sequence_with_errors}\n")
#         outfile.write(f"Number of contigs: {len(contigs)}\n\n")
#         outfile.write("Contigs:\n")
#         for i, contig in enumerate(assembled_contigs, 1):
#             outfile.write(f"Contig {i} (length {len(contig)}):\n")
#             outfile.write(f"{contig}\n")
        
#         # Write longest contig information
#         longest_contig, longest_length = get_longest_contig(assembled_contigs)
#         if longest_contig:
#             outfile.write(f"\nLongest contig (length {longest_length}):\n")
#             outfile.write(f"{longest_contig}\n")
#             outfile.write(f"Percentage of original sequence: {(longest_length/len(sequence))*100:.2f}%\n")
#         else:
#             outfile.write("\nNo contigs found\n")

# if __name__ == "__main__":
#     main()

# from Bio import SeqIO
# import argparse
# import matplotlib.pyplot as plt
# from matplotlib.patches import Circle, FancyArrowPatch
# import math
# from matplotlib.backends.backend_pdf import PdfPages
# from collections import defaultdict, Counter
# import random

# def introduce_errors(sequence, error_rate=0.05, min_gap_size=3):
#     """
#     Introduce errors in the sequence to create gaps/contigs
#     """
#     seq_list = list(sequence)
#     error_positions = []
    
#     # Introduce gaps (removing nucleotides)
#     i = 0
#     while i < len(seq_list):
#         if random.random() < error_rate:
#             gap_size = random.randint(min_gap_size, min_gap_size + 2)
#             if i + gap_size < len(seq_list):
#                 for j in range(gap_size):
#                     seq_list[i + j] = 'N'
#                 error_positions.extend(range(i, i + gap_size))
#                 i += gap_size
#             else:
#                 break
#         i += 1
    
#     return ''.join(seq_list), error_positions

# def find_strongly_connected_components(edges, nodes):
#     """
#     Find strongly connected components in the graph using Kosaraju's algorithm
#     """
#     def dfs_forward(node, visited, order):
#         visited.add(node)
#         if node in edges:
#             for next_node in edges[node]:
#                 if next_node not in visited:
#                     dfs_forward(next_node, visited, order)
#         order.append(node)

#     def dfs_reverse(node, visited, component):
#         visited.add(node)
#         component.append(node)
#         for prev_node in nodes:
#             if prev_node in edges and node in edges[prev_node] and prev_node not in visited:
#                 dfs_reverse(prev_node, visited, component)

#     # First DFS to get order
#     visited = set()
#     order = []
#     for node in nodes:
#         if node not in visited:
#             dfs_forward(node, visited, order)

#     # Second DFS to find SCCs
#     visited = set()
#     components = []
#     for node in reversed(order):
#         if node not in visited:
#             component = []
#             dfs_reverse(node, visited, component)
#             if component:
#                 components.append(component)

#     return components

# def find_contigs(edges, nodes):
#     """
#     Find contigs in the De Bruijn graph by identifying maximal paths
#     """
#     components = find_strongly_connected_components(edges, nodes)
#     contigs = []
    
#     for component in components:
#         # For each component, find maximal non-branching paths
#         visited_edges = set()
        
#         for start_node in component:
#             if start_node not in edges:
#                 continue
                
#             # Check if this node is the start of a non-branching path
#             while start_node in edges and len(edges[start_node]) == 1:
#                 current_path = [start_node]
#                 current = start_node
                
#                 # Follow the path until we reach a branch or visited edge
#                 while current in edges and len(edges[current]) == 1:
#                     next_node = next(iter(edges[current]))
#                     edge = (current, next_node)
                    
#                     if edge in visited_edges:
#                         break
                        
#                     visited_edges.add(edge)
#                     current_path.append(next_node)
#                     current = next_node
                    
#                     # Break if we've returned to start or reached a branch
#                     if current == start_node or (current in edges and len(edges[current]) != 1):
#                         break
                
#                 if len(current_path) > 1:
#                     contigs.append(current_path)
#                 break
    
#     return contigs

# def assemble_contigs(contigs, kmers):
#     """
#     Assemble sequences for each contig
#     """
#     assembled_contigs = []
    
#     for contig in contigs:
#         if len(contig) < 2:  # Skip single-node contigs
#             continue
            
#         # Start with first node
#         sequence = contig[0]
        
#         # Add last character of each subsequent node
#         for node in contig[1:]:
#             sequence += node[-1]
            
#         assembled_contigs.append(sequence)
    
#     return assembled_contigs

# def get_longest_contig(assembled_contigs):
#     """
#     Find the longest contig from the assembled sequences
    
#     Parameters:
#     assembled_contigs (list): List of assembled contig sequences
    
#     Returns:
#     tuple: (longest_contig, length) or (None, 0) if no contigs
#     """
#     if not assembled_contigs:
#         return None, 0
        
#     longest_contig = max(assembled_contigs, key=len)
#     return longest_contig, len(longest_contig)

# def build_debruijn_graph(kmers):
#     """
#     Build a de Bruijn graph from k-mers and their counts
#     """
#     edges = defaultdict(Counter)
#     nodes = set()
    
#     for kmer in kmers:
#         # Skip k-mers containing errors (N's)
#         if 'N' in kmer:
#             continue
            
#         prefix = kmer[:-1]
#         suffix = kmer[1:]
#         nodes.add(prefix)
#         nodes.add(suffix)
#         edges[prefix][suffix] += 1
        
#     return edges, nodes

# def find_kmers(sequence, k):
#     """
#     Find k-mers in sequence order, skipping k-mers containing errors
#     """
#     # Add k-1 characters to handle circular nature
#     extended_sequence = sequence + sequence[:k-1]
    
#     # Get k-mers in order
#     kmers = []
#     for i in range(len(sequence)):
#         kmer = extended_sequence[i:i+k]
#         # Only include k-mers that don't contain errors
#         if 'N' not in kmer:
#             kmers.append(kmer)
    
#     return kmers

# def get_node_positions(nodes):
#     """Calculate node positions in a circle"""
#     positions = {}
#     num_nodes = len(nodes)
#     radius = 1  # Unit circle
    
#     for i, node in enumerate(nodes):
#         angle = 2 * math.pi * i / num_nodes
#         x = radius * math.cos(angle)
#         y = radius * math.sin(angle)
#         positions[node] = (x, y)
    
#     return positions

# def create_figure():
#     """Create and return figure and axis"""
#     fig, ax = plt.subplots(figsize=(12, 12))
#     ax.set_aspect('equal')
#     ax.set_xlim(-1.5, 1.5)
#     ax.set_ylim(-1.5, 1.5)
#     plt.axis('off')
#     return fig, ax

# def draw_nodes(ax, nodes, positions):
#     """Draw nodes as circles with labels"""
#     node_radius = 0.1
#     for node, (x, y) in positions.items():
#         circle = Circle((x, y), node_radius, 
#                       facecolor='lightblue',
#                       edgecolor='black',
#                       alpha=0.6)
#         ax.add_patch(circle)
#         plt.text(x, y, node,
#                 horizontalalignment='center',
#                 verticalalignment='center',
#                 fontweight='bold')
#     return node_radius

# def draw_edges(ax, kmers, positions, node_radius):
#     """Draw edges as arrows with labels"""
#     edge_counts = defaultdict(int)
#     edge_orders = []
    
#     for kmer in kmers:
#         if 'N' not in kmer:  # Only process valid k-mers
#             edge = (kmer[:-1], kmer[1:])
#             edge_counts[edge] += 1
#             edge_orders.append(edge)
    
#     drawn_edges = set()
    
#     for edge in edge_counts:
#         if edge in drawn_edges:
#             continue
            
#         start_node, end_node = edge
#         count = edge_counts[edge]
        
#         if start_node not in positions or end_node not in positions:
#             continue
            
#         start_pos = positions[start_node]
#         end_pos = positions[end_node]
        
#         dx = end_pos[0] - start_pos[0]
#         dy = end_pos[1] - start_pos[1]
#         length = math.sqrt(dx*dx + dy*dy)
        
#         if length == 0:
#             continue

#         start_x = start_pos[0] + (dx * node_radius / length)
#         start_y = start_pos[1] + (dy * node_radius / length)
#         end_x = end_pos[0] - (dx * node_radius / length)
#         end_y = end_pos[1] - (dy * node_radius / length)
        
#         arrow = FancyArrowPatch(
#             (start_x, start_y), (end_x, end_y),
#             arrowstyle='-|>',
#             mutation_scale=20,
#             linewidth=1,
#             color='gray'
#         )
#         ax.add_patch(arrow)
        
#         edge_positions = []
#         for i, e in enumerate(edge_orders):
#             if e == edge:
#                 edge_positions.append(i + 1)
        
#         kmer = start_node + end_node[-1]
#         for i, pos in enumerate(edge_positions):
#             frac = (i + 1)/(count + 1)
#             mid_x = start_x + (end_x - start_x) * frac
#             mid_y = start_y + (end_y - start_y) * frac
            
#             label = f"{pos}:{kmer}"
#             if count > 1:
#                 label += f"(x{count}, #{i+1})"
                
#             plt.text(mid_x, mid_y, label,
#                     bbox=dict(facecolor='white', edgecolor='black', alpha=0.7, pad=2),
#                     horizontalalignment='center',
#                     verticalalignment='center')
        
#         drawn_edges.add(edge)

# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('fasta_file')
#     parser.add_argument('-k', type=int)
#     parser.add_argument('--error_rate', type=float, default=0.0,
#                       help='Rate of error introduction (default: 0.0 - no errors)')
#     args = parser.parse_args()
#     k = args.k
    
#     # Read sequence
#     sequence = str(next(SeqIO.parse(args.fasta_file, "fasta")).seq)
    
#     if len(sequence) < k:
#         print(f"Error: Sequence length ({len(sequence)}) must be >= k ({k})")
#         return
    
#     # Only introduce errors if error_rate is greater than 0
#     if args.error_rate > 0:
#         sequence_with_errors, error_positions = introduce_errors(sequence, args.error_rate)
#         print("\nRunning with introduced errors:")
#         print(f"Error rate: {args.error_rate}")
#         print(f"Number of error positions: {len(error_positions)}")
#     else:
#         sequence_with_errors = sequence
#         error_positions = []
#         print("\nRunning without errors (error_rate = 0)")
    
#     # Find kmers in sequence order
#     kmers = find_kmers(sequence_with_errors, k)
#     kmer_counts = Counter(kmers)
    
#     # Build de Bruijn graph
#     edges, nodes = build_debruijn_graph(kmers)
    
#     # Find contigs
#     contigs = find_contigs(edges, nodes)
#     assembled_contigs = assemble_contigs(contigs, kmers)
    
#     print("\nDe Bruijn Graph Assembly Results:")
#     print(f"k-mer size: {k}")
#     print(f"Original sequence length: {len(sequence)}")
#     print(f"Number of k-mers: {len(kmers)}")
#     print(f"Number of nodes: {len(nodes)}")
#     print(f"Number of contigs: {len(contigs)}")
    
#     if args.error_rate > 0:
#         print(f"\nError information:")
#         print(f"Error rate: {args.error_rate}")
#         print(f"Number of error positions: {len(error_positions)}")
    
#     print("\nContigs:")
#     for i, contig in enumerate(assembled_contigs, 1):
#         print(f"Contig {i} (length {len(contig)}):")
#         print(contig)
    
#     # Get and print longest contig
#     longest_contig, longest_length = get_longest_contig(assembled_contigs)
#     if longest_contig:
#         print(f"\nLongest contig (length {longest_length}):")
#         print(longest_contig)
#     else:
#         print("\nNo contigs found")
    
#     # Calculate node positions and create visualization
#     positions = get_node_positions(list(nodes))
#     fig, ax = create_figure()
#     node_radius = draw_nodes(ax, nodes, positions)
#     draw_edges(ax, kmers, positions, node_radius)
    
#     # Save visualization
#     with PdfPages('debruijn_visual_with_errors.pdf') as pdf:
#         pdf.savefig(fig, bbox_inches='tight')
#     plt.close()
    
#     # Save detailed results to file
#     with open('output_with_errors.txt', 'w') as outfile:
#         outfile.write(f"De Bruijn Graph Assembly Results:\n")
#         outfile.write(f"k-mer size: {k}\n")
#         outfile.write(f"Original sequence length: {len(sequence)}\n")
#         outfile.write(f"Number of k-mers: {len(kmers)}\n")
#         outfile.write(f"Number of nodes: {len(nodes)}\n")
#         outfile.write(f"Original sequence: {sequence}\n")
#         outfile.write(f"Sequence with errors: {sequence_with_errors}\n")
#         outfile.write(f"Number of contigs: {len(contigs)}\n\n")
        
#         if args.error_rate > 0:
#             outfile.write(f"Error information:\n")
#             outfile.write(f"Error rate: {args.error_rate}\n")
#             outfile.write(f"Number of error positions: {len(error_positions)}\n")
#             outfile.write(f"Error positions: {error_positions}\n")
#             outfile.write(f"Original sequence: {sequence}\n")
#             outfile.write(f"Sequence with errors: {sequence_with_errors}\n\n")
#         outfile.write("Contigs:\n")
#         for i, contig in enumerate(assembled_contigs, 1):
#             outfile.write(f"Contig {i} (length {len(contig)}):\n")
#             outfile.write(f"{contig}\n")
        
#         # Write longest contig information
#         longest_contig, longest_length = get_longest_contig(assembled_contigs)
#         if longest_contig:
#             outfile.write(f"\nLongest contig (length {longest_length}):\n")
#             outfile.write(f"{longest_contig}\n")
#             outfile.write(f"Percentage of original sequence: {(longest_length/len(sequence))*100:.2f}%\n")
#         else:
#             outfile.write("\nNo contigs found\n")

# if __name__ == "__main__":
#     main()

"""
This is just for explanation of how the contigs are found and assembled:

-------------------------------------------------------------
link to Kosaraju's Algorithm:
https://en.wikipedia.org/wiki/Kosaraju%27s_algorithm
-------------------------------------------------------------

# Detailed Walkthrough of SCC and Contig Calculation

## Initial Sequence
Original sequence: ATGGCGTGCATGCATGCGTAATC
With errors (20% error rate): ATGGCGTNNNATGCATNNNGTNNTC
k-mer size = 4

## 1. Find Valid K-mers
Skip any k-mers containing 'N':

Valid k-mers (in order):
1. ATGG
2. TGGC
3. GGCG
4. GCGT
5. ATGC
6. TGCA
7. GCAT
8. GTAA
9. TAAC

## 2. Build De Bruijn Graph
Nodes are (k-1)-mers:
- ATG, TGG, GGC, GCG, CGT, TGC, GCA, CAT, GTA, TAA, AAT

Edges (with counts):

edges = {
    'ATG': Counter({'TGG': 1, 'TGC': 1}),
    'TGG': Counter({'GGC': 1}),
    'GGC': Counter({'GCG': 1}),
    'GCG': Counter({'CGT': 1}),
    'TGC': Counter({'GCA': 1}),
    'GCA': Counter({'CAT': 1}),
    'GTA': Counter({'TAA': 1}),
    'TAA': Counter({'AAT': 1})
}

nodes = {'ATG', 'TGG', 'GGC', 'GCG', 'CGT', 'TGC', 'GCA', 'CAT', 'GTA', 'TAA', 'AAT'}

## 3. Finding Strongly Connected Components (SCCs)

### Step-by-Step DFS Forward Pass:

# Initial state
visited = set()
order = []

# Start from node 'ATG'
1. Visit 'ATG'
   - Follow edge to 'TGG'
2. Visit 'TGG'
   - Follow edge to 'GGC'
3. Visit 'GGC'
   - Follow edge to 'GCG'
4. Visit 'GCG'
   - Follow edge to 'CGT'
5. Visit 'CGT'
   - No outgoing edges
   - Add to order: ['CGT']
6. Backtrack to 'GCG'
   - Add to order: ['CGT', 'GCG']
...

Final order: ['CGT', 'GCG', 'GGC', 'TGG', 'GCA', 'CAT', 'TGC', 'TAA', 'AAT', 'GTA', 'ATG']

### Step-by-Step DFS Reverse Pass:
# Process nodes in reverse order
visited = set()
components = []

1. Start with 'ATG':
   Component 1 = ['ATG', 'TGG', 'GGC', 'GCG', 'CGT']
   
2. Start with 'GTA':
   Component 2 = ['GTA', 'TAA', 'AAT']
   
3. Start with 'TGC':
   Component 3 = ['TGC', 'GCA', 'CAT']

Final SCCs = [
    ['ATG', 'TGG', 'GGC', 'GCG', 'CGT'],  # Component 1
    ['GTA', 'TAA', 'AAT'],                 # Component 2
    ['TGC', 'GCA', 'CAT']                  # Component 3
]

## 4. Finding Contigs

### Process Component 1:

visited_edges = set()
current_component = ['ATG', 'TGG', 'GGC', 'GCG', 'CGT']

1. Start at 'ATG':
   - Check edges: len(edges['ATG']) = 2
   - Branching node, skip

2. Start at 'TGG':
   - Check edges: len(edges['TGG']) = 1
   - Follow path:
     TGG -> GGC -> GCG -> CGT
   - No branches, create contig:
     Contig 1 = ['TGGC', 'GGCG', 'GCGT']

3. Continue with remaining nodes...

### Process Component 2:
visited_edges = set()
current_component = ['GTA', 'TAA', 'AAT']

1. Start at 'GTA':
   - Check edges: len(edges['GTA']) = 1
   - Follow path:
     GTA -> TAA -> AAT
   - No branches, create contig:
     Contig 2 = ['GTAA', 'TAAC']

### Process Component 3:
visited_edges = set()
current_component = ['TGC', 'GCA', 'CAT']

1. Start at 'TGC':
   - Check edges: len(edges['TGC']) = 1
   - Follow path:
     TGC -> GCA -> CAT
   - No branches, create contig:
     Contig 3 = ['TGCA', 'GCAT']

## 5. Assemble Contigs
Convert paths to sequences:

Contig 1: TGGCGT   (from ['TGGC', 'GGCG', 'GCGT'])
Contig 2: GTAAC    (from ['GTAA', 'TAAC'])
Contig 3: TGCAT    (from ['TGCA', 'GCAT'])
"""

"""
De Bruijn Graph Assembly Implementation

This script implements a De Bruijn graph-based genome assembly algorithm with error simulation capabilities.
It can process DNA sequences, introduce artificial errors, construct De Bruijn graphs, and identify contigs
through graph traversal.

Key Features:
- Reads DNA sequences from FASTA files
- Optional error introduction to simulate sequencing errors
- De Bruijn graph construction and visualization
- Contig assembly using graph traversal algorithms
- Detailed output generation including assembly statistics
- PDF visualization of the De Bruijn graph structure

Dependencies:
    - BioPython: For reading FASTA files
    - matplotlib: For visualization
    - argparse: For command-line argument parsing
    - collections: For defaultdict and Counter
    - math: For mathematical operations
    - random: For error introduction

Usage:
    python script.py input.fasta -k <kmer_size> [--error_rate <rate>]

Example:
    python script.py sequence.fasta -k 25 --error_rate 0.05
"""

from Bio import SeqIO
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
import math
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict, Counter
import random

def introduce_errors(sequence, error_rate=0.05, min_gap_size=3):
    """
    Introduce artificial errors in a DNA sequence to simulate sequencing errors.
    
    This function creates gaps in the sequence by replacing nucleotides with 'N' characters.
    The gaps are introduced randomly based on the specified error rate and minimum gap size.
    This simulates real-world sequencing errors or low-coverage regions in genome assembly.
    
    Args:
        sequence (str): The original DNA sequence to modify.
        error_rate (float): Probability of introducing an error at each position (default: 0.05).
        min_gap_size (int): Minimum number of consecutive nucleotides to replace (default: 3).
    
    Returns:
        tuple: (
            str: Modified sequence with introduced errors ('N' characters),
            list: Positions where errors were introduced
        )
    
    Example:
        >>> seq = "ATCGATCG"
        >>> modified_seq, error_pos = introduce_errors(seq, error_rate=0.1)
        >>> print(modified_seq)
        'ATCNNNNCG'
    """
    seq_list = list(sequence)
    error_positions = []
    
    i = 0
    while i < len(seq_list):
        if random.random() < error_rate:
            # Calculate gap size: random value between min_gap_size and min_gap_size + 2
            gap_size = random.randint(min_gap_size, min_gap_size + 2)
            
            # Only introduce gap if it fits within remaining sequence
            if i + gap_size < len(seq_list):
                # Replace nucleotides with 'N' characters
                for j in range(gap_size):
                    seq_list[i + j] = 'N'
                # Record positions where errors were introduced
                error_positions.extend(range(i, i + gap_size))
                i += gap_size
            else:
                break
        i += 1
    
    return ''.join(seq_list), error_positions

def dfs_forward(node, edges, visited, order):
    """
    First DFS pass to determine node processing order.
    
    Args:
        node: Current node being processed
        edges (dict): Graph edges {node: set(neighbors)}
        visited (set): Nodes already visited
        order (list): Post-order of node traversal
    """
    visited.add(node)
    if node in edges:
        for next_node in edges[node]:
            if next_node not in visited:
                dfs_forward(next_node, edges, visited, order)
    order.append(node)

def dfs_reverse(node, edges, nodes, visited, component):
    """
    Second DFS pass to identify strongly connected components.
    
    Args:
        node: Current node being processed
        edges (dict): Graph edges {node: set(neighbors)}
        nodes (set): All nodes in the graph
        visited (set): Nodes already visited
        component (list): Current component being built
    """
    visited.add(node)
    component.append(node)
    for prev_node in nodes:
        if prev_node in edges and node in edges[prev_node] and prev_node not in visited:
            dfs_reverse(prev_node, edges, nodes, visited, component)

def find_strongly_connected_components(edges, nodes):
    """
    Identify strongly connected components (SCCs) in the De Bruijn graph using Kosaraju's algorithm.
    
    This implementation uses a two-pass depth-first search approach:
    1. First pass: DFS to get finishing times (reverse post-order)
    2. Second pass: DFS on transpose graph to find SCCs
    
    SCCs are important in genome assembly as they can represent:
    - Repeated regions in the genome
    - Complex graph structures that need special handling
    - Potential assembly challenges
    
    Args:
        edges (dict): Adjacency list representation of graph edges {node: set(neighbors)}.
        nodes (set): Set of all nodes in the graph.
    
    Returns:
        list: Lists of nodes forming strongly connected components.
    """
    # First DFS pass: get node processing order
    visited = set()
    order = []
    for node in nodes:
        if node not in visited:
            dfs_forward(node, edges, visited, order)

    # Second DFS pass: find SCCs
    visited = set()
    components = []
    # Process nodes in reverse order of finishing times
    for node in reversed(order):
        if node not in visited:
            component = []
            dfs_reverse(node, edges, nodes, visited, component)
            if component:  # Only add non-empty components
                components.append(component)

    return components

def find_contigs(edges, nodes):
    """
    Identify contigs in the De Bruijn graph by finding maximal non-branching paths.
    
    A contig is a contiguous sequence of DNA that can be reliably assembled. In the
    context of a De Bruijn graph, contigs are represented by maximal non-branching
    paths - paths that don't contain any nodes with multiple incoming or outgoing edges
    (except possibly at the ends).
    
    The algorithm:
    1. Finds strongly connected components
    2. Within each component, identifies maximal non-branching paths
    3. Tracks visited edges to avoid duplicate paths
    
    Args:
        edges (dict): Adjacency list representation of graph edges {node: set(neighbors)}.
        nodes (set): Set of all nodes in the graph.
    
    Returns:
        list: Lists of nodes representing contigs (maximal non-branching paths).
    
    Note:
        - A path is non-branching if each internal node has exactly one incoming and one outgoing edge
        - The algorithm handles both linear and circular paths in the graph
    """
    components = find_strongly_connected_components(edges, nodes)
    contigs = []
    
    for component in components:
        visited_edges = set()
        
        for start_node in component:
            if start_node not in edges:
                continue
                
            # Check if this node could be the start of a non-branching path
            while start_node in edges and len(edges[start_node]) == 1:
                current_path = [start_node]
                current = start_node
                
                # Follow the path until we reach a branch or visited edge
                while current in edges and len(edges[current]) == 1:
                    next_node = next(iter(edges[current]))
                    edge = (current, next_node)
                    
                    # Stop if we've seen this edge before
                    if edge in visited_edges:
                        break
                        
                    visited_edges.add(edge)
                    current_path.append(next_node)
                    current = next_node
                    
                    # Stop if we've completed a cycle or reached a branching point
                    if current == start_node or (current in edges and len(edges[current]) != 1):
                        break
                
                # Only keep paths with at least 2 nodes
                if len(current_path) > 1:
                    contigs.append(current_path)
                break
    
    return contigs

def assemble_contigs(contigs, kmers):
    """
    Convert contig paths into DNA sequences.
    
    This function takes the paths found in the De Bruijn graph (contigs) and
    reconstructs the corresponding DNA sequences. It uses the overlap property
    of k-mers to efficiently join the sequences.
    
    Assembly process:
    1. Start with the first node's sequence
    2. For each subsequent node, append only its last character
    3. This works because adjacent nodes overlap by k-1 characters
    
    Args:
        contigs (list): Lists of nodes representing contig paths.
        kmers (list): Original k-mers used to build the graph.
    
    Returns:
        list: DNA sequences corresponding to the assembled contigs.
    
    Example:
        If k=3 and a contig path is ['ATG', 'TGC', 'GCT']:
        The assembled sequence would be 'ATGCT'
        (because 'ATG' overlaps with 'TGC' by 'TG', and 'TGC' overlaps with 'GCT' by 'GC')
    """
    assembled_contigs = []
    
    for contig in contigs:
        if len(contig) < 2:  # Skip single-node contigs
            continue
            
        # Start with the complete sequence of the first node
        sequence = contig[0]
        
        # For each subsequent node, we only need to append its last character
        # since it overlaps with the previous node by k-1 characters
        for node in contig[1:]:
            sequence += node[-1]
            
        assembled_contigs.append(sequence)
    
    return assembled_contigs

def get_longest_contig(assembled_contigs):
    """
    Find the longest contig among the assembled sequences.
    
    In genome assembly, longer contigs are generally more reliable and
    provide more useful information. This function identifies the longest
    contig and its length from the assembly results.
    
    Args:
        assembled_contigs (list): List of assembled contig sequences.
    
    Returns:
        tuple: (
            str or None: The longest contig sequence (None if no contigs),
            int: Length of the longest contig (0 if no contigs)
        )
    """
    if not assembled_contigs:
        return None, 0
        
    longest_contig = max(assembled_contigs, key=len)
    return longest_contig, len(longest_contig)

def build_debruijn_graph(kmers):
    """
    Construct a De Bruijn graph from a set of k-mers.
    
    A De Bruijn graph is a directed graph where:
    - Nodes are (k-1)-mers
    - Edges represent k-mers, connecting overlapping (k-1)-mers
    - Edge weights (counts) track k-mer frequencies
    
    Construction process:
    1. For each k-mer, create nodes for prefix and suffix (k-1)-mers
    2. Add directed edge from prefix to suffix
    3. Track edge multiplicities using Counter
    
    Args:
        kmers (list): List of k-mers from the input sequence.
    
    Returns:
        tuple: (
            dict: Edge dictionary {prefix: Counter(suffixes)},
            set: Set of all nodes (unique (k-1)-mers)
        )
    
    Note:
        - Skips k-mers containing errors (marked by 'N')
        - Edge counts are important for:
            * Handling repeats in the genome
            * Identifying potential sequencing errors
            * Resolving assembly ambiguities
    """
    edges = defaultdict(Counter)
    nodes = set()
    
    for kmer in kmers:
        # Skip k-mers containing errors
        if 'N' in kmer:
            continue
            
        # Each k-mer becomes an edge connecting its prefix and suffix
        prefix = kmer[:-1]  # First k-1 characters
        suffix = kmer[1:]   # Last k-1 characters
        
        # Add nodes and edge
        nodes.add(prefix)
        nodes.add(suffix)
        edges[prefix][suffix] += 1
        
    return edges, nodes

def find_kmers(sequence, k):
    """
    Extract all valid k-mers from a DNA sequence.
    
    This function handles circular sequences by extending the original sequence
    with k-1 characters from the start. It also filters out k-mers containing
    error characters ('N').
    
    Args:
        sequence (str): Input DNA sequence.
        k (int): Size of k-mers to extract.
    
    Returns:
        list: Valid k-mers in their order of appearance in the sequence.
    
    Note:
        - Handles circular sequences by wrapping around
        - Skips k-mers containing 'N' (error positions)
        - Maintains original sequence order, important for:
            * Tracking k-mer positions
            * Preserving sequence context
            * Assembly validation
    """
    # Handle circular nature of sequence by adding k-1 chars from start
    extended_sequence = sequence + sequence[:k-1]
    
    # Extract k-mers, skipping those with errors
    kmers = []
    for i in range(len(sequence)):
        kmer = extended_sequence[i:i+k]
        if 'N' not in kmer:
            kmers.append(kmer)
    
    return kmers

def get_node_positions(nodes):
    """
    Calculate node positions for circular graph visualization.
    
    Positions nodes evenly around a unit circle for aesthetic graph drawing.
    This layout works well for De Bruijn graphs as it:
    - Minimizes edge crossings
    - Provides clear visualization of graph structure
    - Scales well with different numbers of nodes
    
    Args:
        nodes (set): Set of all nodes to position.
    
    Returns:
        dict: Mapping of nodes to (x,y) coordinates.
    """
    positions = {}
    num_nodes = len(nodes)
    radius = 1  # Unit circle
    
    for i, node in enumerate(nodes):
        angle = 2 * math.pi * i / num_nodes
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        positions[node] = (x, y)
    
    return positions

def create_figure():
    """
    Initialize matplotlib figure for graph visualization.
    
    Creates a square figure with appropriate dimensions and settings for
    the De Bruijn graph visualization. The visualization is centered on
    the origin with equal scaling on both axes.
    
    Returns:
        tuple: (matplotlib.figure.Figure, matplotlib.axes.Axes)
    """
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.set_aspect('equal')
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    plt.axis('off')
    return fig, ax

def draw_nodes(ax, nodes, positions):
    """
    Draw graph nodes as circles with labels.
    
    Creates a visual representation of each node in the De Bruijn graph:
    - Nodes are drawn as circles with a light blue fill
    - Each node is labeled with its sequence
    - Consistent sizing ensures readability
    
    Args:
        ax (matplotlib.axes.Axes): The figure axis for drawing
        nodes (set): Set of all nodes to draw
        positions (dict): Mapping of nodes to (x,y) coordinates
    
    Returns:
        float: Radius used for the nodes (needed for edge drawing)
    """
    node_radius = 0.1
    for node, (x, y) in positions.items():
        # Create and add circle for node
        circle = Circle((x, y), node_radius, 
                      facecolor='lightblue',
                      edgecolor='black',
                      alpha=0.6)
        ax.add_patch(circle)
        
        # Add node label (sequence)
        plt.text(x, y, node,
                horizontalalignment='center',
                verticalalignment='center',
                fontweight='bold')
    return node_radius

def draw_edges(ax, kmers, positions, node_radius):
    """
    Draw graph edges as arrows with detailed labels.
    
    Creates a visual representation of the edges in the De Bruijn graph:
    - Edges are drawn as arrows connecting nodes
    - Labels show:
        * Edge order in the original sequence
        * The complete k-mer sequence
        * Multiplicity information for repeated edges
    
    Args:
        ax (matplotlib.axes.Axes): The figure axis for drawing
        kmers (list): List of k-mers representing edges
        positions (dict): Mapping of nodes to (x,y) coordinates
        node_radius (float): Radius of nodes (for proper arrow positioning)
    
    Note:
        Edge labels include:
        - Position number in sequence
        - k-mer sequence
        - Count information for repeated edges
        - Multiple labels for edges appearing multiple times
    """
    edge_counts = defaultdict(int)
    edge_orders = []
    
    # First pass: Count edges and record their order
    for kmer in kmers:
        if 'N' not in kmer:  # Skip error-containing k-mers
            edge = (kmer[:-1], kmer[1:])
            edge_counts[edge] += 1
            edge_orders.append(edge)
    
    drawn_edges = set()
    
    # Second pass: Draw edges and labels
    for edge in edge_counts:
        if edge in drawn_edges:
            continue
            
        start_node, end_node = edge
        count = edge_counts[edge]
        
        # Skip edges with missing positions
        if start_node not in positions or end_node not in positions:
            continue
            
        start_pos = positions[start_node]
        end_pos = positions[end_node]
        
        # Calculate edge vector
        dx = end_pos[0] - start_pos[0]
        dy = end_pos[1] - start_pos[1]
        length = math.sqrt(dx*dx + dy*dy)
        
        # Skip self-loops
        if length == 0:
            continue

        # Adjust arrow endpoints to account for node sizes
        start_x = start_pos[0] + (dx * node_radius / length)
        start_y = start_pos[1] + (dy * node_radius / length)
        end_x = end_pos[0] - (dx * node_radius / length)
        end_y = end_pos[1] - (dy * node_radius / length)
        
        # Draw the arrow
        arrow = FancyArrowPatch(
            (start_x, start_y), (end_x, end_y),
            arrowstyle='-|>',
            mutation_scale=20,
            linewidth=1,
            color='gray'
        )
        ax.add_patch(arrow)
        
        # Find all positions where this edge appears
        edge_positions = []
        for i, e in enumerate(edge_orders):
            if e == edge:
                edge_positions.append(i + 1)
        
        # Create and place labels for each occurrence
        kmer = start_node + end_node[-1]
        for i, pos in enumerate(edge_positions):
            # Calculate label position along edge
            frac = (i + 1)/(count + 1)
            mid_x = start_x + (end_x - start_x) * frac
            mid_y = start_y + (end_y - start_y) * frac
            
            # Construct label text
            label = f"{pos}:{kmer}"
            if count > 1:
                label += f"(x{count}, #{i+1})"
                
            # Draw label with background box
            plt.text(mid_x, mid_y, label,
                    bbox=dict(facecolor='white', edgecolor='black', alpha=0.7, pad=2),
                    horizontalalignment='center',
                    verticalalignment='center')
        
        drawn_edges.add(edge)

def main():
    """
    Main function handling the complete De Bruijn graph assembly workflow.
    
    Workflow steps:
    1. Parse command line arguments
    2. Read and validate input sequence
    3. Optionally introduce errors
    4. Build and analyze De Bruijn graph
    5. Find and assemble contigs
    6. Generate visualization
    7. Save detailed results
    
    Command-line arguments:
    - fasta_file: Input FASTA file containing sequence
    - k: k-mer size for graph construction
    - error_rate: Optional rate for introducing errors
    
    Outputs:
    - Console summary of assembly results
    - PDF visualization of the De Bruijn graph
    - Detailed text output file with assembly statistics
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file',
                       help='Input FASTA file containing sequence to assemble')
    parser.add_argument('-k', type=int,
                       help='Size of k-mers for De Bruijn graph construction')
    parser.add_argument('--error_rate', type=float, default=0.0,
                       help='Rate of error introduction (default: 0.0 - no errors)')
    args = parser.parse_args()
    k = args.k
    
    # Read input sequence
    sequence = str(next(SeqIO.parse(args.fasta_file, "fasta")).seq)
    
    # Validate sequence length
    if len(sequence) < k:
        print(f"Error: Sequence length ({len(sequence)}) must be >= k ({k})")
        return
    
    # Handle error introduction
    if args.error_rate > 0:
        sequence_with_errors, error_positions = introduce_errors(sequence, args.error_rate)
        print("\nRunning with introduced errors:")
        print(f"Error rate: {args.error_rate}")
        print(f"Number of error positions: {len(error_positions)}")
    else:
        sequence_with_errors = sequence
        error_positions = []
        print("\nRunning without errors (error_rate = 0)")
    
    # Extract k-mers and build graph
    kmers = find_kmers(sequence_with_errors, k)
    kmer_counts = Counter(kmers)
    edges, nodes = build_debruijn_graph(kmers)
    
    # Find and assemble contigs
    contigs = find_contigs(edges, nodes)
    assembled_contigs = assemble_contigs(contigs, kmers)
    
    # Print assembly statistics
    print("\nDe Bruijn Graph Assembly Results:")
    print(f"k-mer size: {k}")
    print(f"Original sequence length: {len(sequence)}")
    print(f"Number of k-mers: {len(kmers)}")
    print(f"Number of nodes: {len(nodes)}")
    print(f"Number of contigs: {len(contigs)}")
    
    # Print error statistics if applicable
    if args.error_rate > 0:
        print(f"\nError information:")
        print(f"Error rate: {args.error_rate}")
        print(f"Number of error positions: {len(error_positions)}")
    
    # Print contig information
    print("\nContigs:")
    for i, contig in enumerate(assembled_contigs, 1):
        print(f"Contig {i} (length {len(contig)}):")
        print(contig)
    
    # Print longest contig information
    longest_contig, longest_length = get_longest_contig(assembled_contigs)
    if longest_contig:
        print(f"\nLongest contig (length {longest_length}):")
        print(longest_contig)
    else:
        print("\nNo contigs found because no valid k-mers foudn for specified k-mer size (error rate too high)")
    
    # Generate visualization
    positions = get_node_positions(list(nodes))
    fig, ax = create_figure()
    node_radius = draw_nodes(ax, nodes, positions)
    draw_edges(ax, kmers, positions, node_radius)
    
    # Save visualization to PDF
    with PdfPages('debruijn_visual_with_errors.pdf') as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    
    # Save detailed results to text file
    with open('output_with_errors.txt', 'w') as outfile:
        outfile.write(f"De Bruijn Graph Assembly Results:\n")
        outfile.write(f"k-mer size: {k}\n")
        outfile.write(f"Original sequence length: {len(sequence)}\n")
        outfile.write(f"Number of k-mers: {len(kmers)}\n")
        outfile.write(f"No. of unique k-mers: {len(kmer_counts)}\n")
        outfile.write(f"Number of nodes: {len(nodes)}\n")
        outfile.write(f"Original sequence: {sequence}\n")
        outfile.write(f"Sequence with errors: {sequence_with_errors}\n")
        outfile.write(f"Same Sequence: {sequence == sequence_with_errors}\n\n")
        outfile.write(f"Number of contigs: {len(contigs)}\n\n")
        
        # Write error information if applicable
        if args.error_rate > 0:
            outfile.write(f"Error information:\n")
            outfile.write(f"Error rate: {args.error_rate}\n")
            outfile.write(f"Number of error positions: {len(error_positions)}\n")
            outfile.write(f"Error positions: {error_positions}\n")
            outfile.write(f"Original sequence: {sequence}\n")
            outfile.write(f"Assembled Sequence: {sequence_with_errors}\n\n")

        # Write contig information
        outfile.write("Contigs:\n")
        for i, contig in enumerate(assembled_contigs, 1):
            outfile.write(f"Contig {i} (length {len(contig)}):\n")
            outfile.write(f"{contig}\n")
        
        # Write longest contig information
        longest_contig, longest_length = get_longest_contig(assembled_contigs)
        if longest_contig:
            outfile.write(f"\nLongest contig (length {longest_length}):\n")
            outfile.write(f"{longest_contig}\n")
            outfile.write(f"Percentage of original sequence: {(longest_length/len(sequence))*100:.2f}%\n")
        else:
            outfile.write("\nNo contigs found because no valid k-mers foudn for specified k-mer size (error rate too high)\n")

if __name__ == "__main__":
    main()