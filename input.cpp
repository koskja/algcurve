#include "input.hpp"
#include <stdexcept>
#include <cctype>

std::string make_error_message(std::string_view expression, size_t pos, const std::string& message) {
    const size_t context_size = 10;  // How many characters of context to show on EACH side
    size_t start = (pos > context_size) ? pos - context_size : 0;
    size_t end = std::min(expression.length(), pos + context_size);
    
    std::string result = message + "\n";
    result += "Position " + std::to_string(pos) + ":\n";
    
    // Add the context with the error position marked
    result += std::string(expression.substr(start, end - start)) + "\n";
    result += std::string(pos - start, ' ') + "^";
    
    return result;
}

std::vector<Token> tokenize(std::string_view expression) {
    std::vector<Token> tokens;
    size_t pos = 0;

    while (pos < expression.length()) {
        while (pos < expression.length() && std::isspace(expression[pos])) {
            pos++;
        }
        if (pos >= expression.length()) break;

        char c = expression[pos];
        TokenType type;
        int onechar = 1;
        switch (c) {
            case '+': type = TokenType::Add; break;
            case '-': type = TokenType::Sub; break;
            case '*': type = TokenType::Mul; break;
            case '^': type = TokenType::Pow; break;
            case '(': type = TokenType::LPar; break;
            case ')': type = TokenType::RPar; break;
            default: onechar = 0; break;
        }
        if (onechar) {
            tokens.push_back({type, expression.substr(pos, 1)});
            pos++;
        }
        else if (std::isalpha(c)) {
            tokens.push_back({TokenType::Var, expression.substr(pos, 1)});
            pos++;
        }
        else if (std::isdigit(c) || c == '.') {
            size_t start = pos++;
            bool hasDecimal = (c == '.');
            
            while (pos < expression.length()) {
                if (std::isdigit(expression[pos])) {
                    pos++;
                }
                else if (expression[pos] == '.' && !hasDecimal) {
                    hasDecimal = true;
                    pos++;
                }
                else break;
            }
            
            if (hasDecimal && pos - start == 1) {
                throw std::runtime_error(make_error_message(expression, start,
                    "Invalid number format: lone decimal point"));
            }
            
            tokens.push_back({TokenType::Const, expression.substr(start, pos - start)});
        }
        else {
            throw std::runtime_error(make_error_message(expression, pos,
                "Invalid character encountered: " + std::string(1, c)));
        }
    }
    
    return tokens;
}

std::unique_ptr<ASTNode> parse_tokens(std::vector<Token> tokens) {
    std::vector<std::unique_ptr<ASTNode>> output_stack;
    std::vector<TokenType> operator_stack;

    auto make_node = [](TokenType op, std::unique_ptr<ASTNode> right, std::unique_ptr<ASTNode> left) {
        auto node = std::make_unique<ASTNode>();
        switch (op) {
            case TokenType::Add:
                node->type = ASTNodeType::Add;
                break;
            case TokenType::Sub:
                node->type = ASTNodeType::Sub;
                break;
            case TokenType::Mul:
                node->type = ASTNodeType::Mul;
                break;
            case TokenType::Pow:
                node->type = ASTNodeType::Pow;
                break;
            default:
                throw std::runtime_error("Invalid operator");
        }
        
        node->value = ASTNodeForkType{std::move(left), std::move(right)};
        return node;
    };

    auto precedence = [](TokenType op) -> int {
        switch (op) {
            case TokenType::Add: return 1;
            case TokenType::Sub: return 1;
            case TokenType::Mul: return 2;
            case TokenType::Pow: return 3;
            default: return 0;
        }
    };

    auto process_operator = [&]() {
        auto op = operator_stack.back();
        operator_stack.pop_back();
        
        if (output_stack.size() < 2) {
            throw std::runtime_error("Invalid expression: not enough operands");
        }
        
        auto right = std::move(output_stack.back());
        output_stack.pop_back();
        auto left = std::move(output_stack.back());
        output_stack.pop_back();
        
        output_stack.push_back(make_node(op, std::move(right), std::move(left)));
    };

    for (size_t i = 0; i < tokens.size(); i++) {
        const auto& token = tokens[i];
        
        switch (token.type) {
            case TokenType::Const: {
                auto node = std::make_unique<ASTNode>();
                node->type = ASTNodeType::Const;
                node->value = std::stod(std::string(token.value));
                output_stack.push_back(std::move(node));
                
                // Handle implicit multiplication with next token
                if (i + 1 < tokens.size() && 
                    (tokens[i + 1].type == TokenType::Var || 
                     tokens[i + 1].type == TokenType::LPar || 
                     tokens[i + 1].type == TokenType::Const)) {
                    while (!operator_stack.empty() && precedence(operator_stack.back()) >= precedence(TokenType::Mul)) {
                        process_operator();
                    }
                    operator_stack.push_back(TokenType::Mul);
                }
                break;
            }
            case TokenType::Var: {
                auto node = std::make_unique<ASTNode>();
                node->type = ASTNodeType::Var;
                node->value = token.value[0];
                output_stack.push_back(std::move(node));
                
                // Handle implicit multiplication with next token
                if (i + 1 < tokens.size() && 
                    (tokens[i + 1].type == TokenType::Var || 
                     tokens[i + 1].type == TokenType::LPar || 
                     tokens[i + 1].type == TokenType::Const)) {
                    while (!operator_stack.empty() && precedence(operator_stack.back()) >= precedence(TokenType::Mul)) {
                        process_operator();
                    }
                    operator_stack.push_back(TokenType::Mul);
                }
                break;
            }
            case TokenType::Add:
            case TokenType::Sub:
            case TokenType::Mul: {
                while (!operator_stack.empty() && precedence(operator_stack.back()) >= precedence(token.type)) {
                    process_operator();
                }
                operator_stack.push_back(token.type);
                break;
            }
            case TokenType::Pow: {
                while (!operator_stack.empty() && precedence(operator_stack.back()) > precedence(token.type)) {
                    process_operator();
                }
                operator_stack.push_back(token.type);
                break;
            }
            case TokenType::LPar:
                operator_stack.push_back(token.type);
                break;
            case TokenType::RPar: {
                while (!operator_stack.empty() && operator_stack.back() != TokenType::LPar) {
                    process_operator();
                }
                if (operator_stack.empty()) {
                    throw std::runtime_error("Mismatched parentheses");
                }
                operator_stack.pop_back(); // Remove LPar
                
                // Handle implicit multiplication after closing parenthesis
                if (i + 1 < tokens.size() && 
                    (tokens[i + 1].type == TokenType::Var || 
                     tokens[i + 1].type == TokenType::LPar || 
                     tokens[i + 1].type == TokenType::Const)) {
                    while (!operator_stack.empty() && precedence(operator_stack.back()) >= precedence(TokenType::Mul)) {
                        process_operator();
                    }
                    operator_stack.push_back(TokenType::Mul);
                }
                break;
            }
        }
    }

    while (!operator_stack.empty()) {
        if (operator_stack.back() == TokenType::LPar) {
            throw std::runtime_error("Mismatched parentheses");
        }
        process_operator();
    }

    if (output_stack.empty()) {
        throw std::runtime_error("Empty expression");
    }
    if (output_stack.size() > 1) {
        throw std::runtime_error("Invalid expression: too many operands");
    }

    return std::move(output_stack.back());
}