
import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from difflib import SequenceMatcher
from io import StringIO

st.set_page_config(
    page_title="DNA + Protein Mutation Analyzer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 DNA + Protein Mutation Analyzer")
st.write(
    "Compare two coding DNA sequences and automatically "
    "analyze their nucleotide and protein differences."
)

st.info(
    "Educational tool: enter sequences in the same coding "
    "orientation and reading frame. Results are not clinical diagnoses."
)

def clean_dna(sequence):
    return "".join(sequence.upper().split())

def valid_dna(sequence):
    return bool(sequence) and all(
        base in "ATGC" for base in sequence
    )

def translate_dna(dna):
    # Translate complete codons only
    usable = dna[:len(dna) // 3 * 3]
    return str(Seq(usable).translate())

def dna_stats(dna):
    length = len(dna)
    counts = {base: dna.count(base) for base in "ATGC"}
    gc = (counts["G"] + counts["C"]) / length * 100
    return length, counts, gc

def dna_changes(ref, var):
    matcher = SequenceMatcher(
        None, ref, var, autojunk=False
    )
    rows = []

    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            continue

        if tag == "replace":
            kind = "Substitution/replacement"
            position = f"{i1 + 1}-{i2}"
            old = ref[i1:i2]
            new = var[j1:j2]
        elif tag == "insert":
            kind = "Insertion"
            position = f"After position {i1}"
            old = "-"
            new = var[j1:j2]
        else:
            kind = "Deletion"
            position = f"{i1 + 1}-{i2}"
            old = ref[i1:i2]
            new = "-"

        rows.append({
            "Change": kind,
            "Reference position": position,
            "Reference bases": old,
            "Variant bases": new
        })

    return rows

def classify_codon(old, new):
    old_aa = str(Seq(old).translate())
    new_aa = str(Seq(new).translate())

    if old_aa == new_aa:
        return "Synonymous (no amino acid change)"
    if old_aa == "*" and new_aa != "*":
        return "Stop-loss"
    if old_aa != "*" and new_aa == "*":
        return "Stop-gain (nonsense)"
    return "Missense (amino acid change)"

def protein_changes(ref_protein, var_protein):
    matcher = SequenceMatcher(
        None, ref_protein, var_protein, autojunk=False
    )
    rows = []

    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            continue

        rows.append({
            "Change": tag.capitalize(),
            "Reference amino acid positions": (
                f"{i1 + 1}-{i2}" if i1 < i2 else f"After {i1}"
            ),
            "Reference amino acids": ref_protein[i1:i2] or "-",
            "Variant amino acids": var_protein[j1:j2] or "-"
        })

    return rows

def codon_analysis(ref, var):
    results = []

    if len(ref) != len(var):
        return results

    usable = len(ref) // 3 * 3

    for i in range(0, usable, 3):
        old = ref[i:i+3]
        new = var[i:i+3]

        if old != new:
            results.append({
                "Codon position": i // 3 + 1,
                "DNA position": f"{i+1}-{i+3}",
                "Reference codon": old,
                "Variant codon": new,
                "Reference amino acid": str(Seq(old).translate()),
                "Variant amino acid": str(Seq(new).translate()),
                "Classification": classify_codon(old, new)
            })

    return results

with st.sidebar:
    st.header("How to use")
    st.write("1. Paste reference DNA.")
    st.write("2. Paste variant DNA.")
    st.write("3. Click Analyze sequences.")
    st.write("4. Download the report.")
    st.caption("Use A, T, G and C only.")

col1, col2 = st.columns(2)

with col1:
    reference_input = st.text_area(
        "Reference DNA sequence",
        value="ATGGCTGAATTTTAA",
        height=150,
        help="Your original coding DNA sequence."
    )

with col2:
    variant_input = st.text_area(
        "Variant DNA sequence",
        value="ATGGTTGAATTTTAA",
        height=150,
        help="The DNA sequence you want to compare."
    )

if st.button("🔬 Analyze sequences", type="primary"):
    ref = clean_dna(reference_input)
    var = clean_dna(variant_input)

    if not valid_dna(ref) or not valid_dna(var):
        st.error("Invalid input. Use DNA sequences containing only A, T, G and C.")
    else:
        ref_protein = translate_dna(ref)
        var_protein = translate_dna(var)

        st.header("1. DNA analysis")

        a, b, c = st.columns(3)
        a.metric("Reference length", f"{len(ref)} nt")
        b.metric("Variant length", f"{len(var)} nt")
        c.metric("Length difference", f"{len(var) - len(ref):+} nt")

        ref_len, ref_counts, ref_gc = dna_stats(ref)
        var_len, var_counts, var_gc = dna_stats(var)

        stats = pd.DataFrame(
            [ref_counts, var_counts],
            index=["Reference", "Variant"]
        )
        stats["Length"] = [ref_len, var_len]
        stats["GC content (%)"] = [round(ref_gc, 2), round(var_gc, 2)]

        st.dataframe(stats, use_container_width=True)

        st.header("2. DNA differences")
        dna_rows = dna_changes(ref, var)

        if dna_rows:
            st.dataframe(pd.DataFrame(dna_rows), use_container_width=True)
        else:
            st.success("The DNA sequences are identical.")

        st.header("3. Protein translation")

        if len(ref) % 3 or len(var) % 3:
            st.warning(
                "At least one sequence has an incomplete codon. "
                "The final 1–2 bases are excluded from translation."
            )

        st.write("**Reference protein**")
        st.code(ref_protein or "(no complete codons)")

        st.write("**Variant protein**")
        st.code(var_protein or "(no complete codons)")

        st.header("4. Protein differences")
        protein_rows = protein_changes(ref_protein, var_protein)

        if protein_rows:
            st.dataframe(
                pd.DataFrame(protein_rows),
                use_container_width=True
            )
        else:
            st.success("The translated protein sequences are identical.")

        st.header("5. Automatic mutation interpretation")

        if len(ref) != len(var):
            difference = len(var) - len(ref)

            if difference % 3 != 0:
                st.warning(
                    "Possible frameshift: sequence lengths differ by "
                    "a number of nucleotides not divisible by three. "
                    "Confirm using a biological sequence alignment."
                )
            else:
                st.info(
                    "Length difference is divisible by three. "
                    "This is compatible with an in-frame insertion or "
                    "deletion, but the exact event needs alignment."
                )
        else:
            codon_rows = codon_analysis(ref, var)

            if not codon_rows:
                st.success("No codon differences detected.")
            else:
                st.dataframe(
                    pd.DataFrame(codon_rows),
                    use_container_width=True
                )

                types = sorted({
                    row["Classification"] for row in codon_rows
                })
                st.write("**Detected codon-level outcomes:**")
                for item in types:
                    st.write(f"- {item}")

        if ref_protein == var_protein:
            st.info(
                "Overall result: no difference in the translated "
                "amino acid sequence."
            )
        else:
            st.write(
                "**Overall result:** the translated protein sequences differ."
            )

        if "*" in var_protein and "*" not in ref_protein:
            st.warning(
                "A stop symbol appears in the variant protein but not "
                "the reference. Check whether a stop-gain or other "
                "sequence change explains this result."
            )

        st.caption(
            "Interpretation assumes coding-strand input, the correct "
            "reading frame, and the standard genetic code. Complex "
            "insertions/deletions and repeat regions may need proper "
            "sequence alignment. A sequence difference alone does "
            "not establish disease causation."
        )

        report = StringIO()
        report.write("DNA + PROTEIN MUTATION ANALYZER REPORT\n\n")
        report.write(f"Reference DNA: {ref}\n")
        report.write(f"Variant DNA: {var}\n")
        report.write(f"Reference protein: {ref_protein}\n")
        report.write(f"Variant protein: {var_protein}\n")
        report.write(f"Reference GC content: {ref_gc:.2f}%\n")
        report.write(f"Variant GC content: {var_gc:.2f}%\n\n")
        report.write("DNA differences:\n")
        report.write(pd.DataFrame(dna_rows).to_string(index=False)
                     if dna_rows else "None")
        report.write("\n\nProtein differences:\n")
        report.write(pd.DataFrame(protein_rows).to_string(index=False)
                     if protein_rows else "None")
        report.write("\n")

        st.download_button(
            "📄 Download results report",
            data=report.getvalue(),
            file_name="mutation_analysis_report.txt",
            mime="text/plain"
        )
