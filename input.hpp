#pragma once
#include <memory>
#include <variant>
#include <vector>
#include <string_view>

enum class TokenType {
    Add,
    Sub,
    Mul,
    Pow,
    Var,
    Const,
    LPar,
    RPar,
};

struct Token {
    TokenType type;
    std::string_view value;
};

enum class ASTNodeType {
    Add,
    Sub,
    Mul,
    Pow,
    Var,
    Const,
};

struct ASTNode;
using ASTNodePtrType = std::unique_ptr<ASTNode>;
using ASTNodeForkType = std::pair<ASTNodePtrType, ASTNodePtrType>;
using ASTNodeConstType = double;
using ASTNodeVarType = char;

struct ASTNode {
    ASTNodeType type;
    std::variant<ASTNodeForkType, ASTNodeVarType, ASTNodeConstType> value;
};

std::unique_ptr<ASTNode> parse_tokens(std::vector<Token> tokens);
std::vector<Token> tokenize(std::string_view expression);
